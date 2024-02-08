"""
Helper functions for preforming/displaying static simulations of the CBP butanol model.
"""

from cobra.flux_analysis import pfba, flux_variability_analysis, production_envelope
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_productions(model, medium: dict = None, reactions: list = None):
    """Get the pFBA solution and (loopless) FVA solution for specified reactions in a model."""
    
    if reactions is None:
        reactions = ["EX_but_e", "EX_ac_e", "EX_etoh_e", "EX_btoh_e", "EX_acetone_e"]

    with model:

        if medium is not None:
            model.medium = medium

        pfba_sol = pfba(model)
        fva_sol = flux_variability_analysis(model, reactions, loopless=True)

        output_df = pd.concat([pfba_sol.fluxes, fva_sol], axis=1, join="inner")
        return output_df.rename(columns={"fluxes": "pFBA sol"})
    

def calculate_essential_rx(model):
    """Calculate essential reactions and metabolites in a model."""

    medium = list(model.medium.keys())
    essential_rx = []

    for rx in medium:
        with model:
            model.reactions.get_by_id(rx).lower_bound = 0
            sol = model.slim_optimize()
            if sol < 1e-6:
                essential_rx.append(rx)

    print('essential metabolites:' ,[rx[3:] for rx in essential_rx])
    print('essential reactions:' ,[rx for rx in essential_rx])


def read_production_rates(sol, production_reactions = None):
    """Returns some metrics from a cobra solution in an easily readable format."""
    if production_reactions is None:
        production_reactions = ["EX_but_e", "EX_ac_e", "EX_etoh_e", "EX_btoh_e", "EX_acetone_e"]
    print('Production fluxes:')
    print(sol.fluxes[production_reactions])


def get_specific_medium(model, specific_reactions: dict, fill_value: float = 0.1):
    """Get a dictionary of max uptake rates for a model, with a set value (default 0.1) for each reaction in the original medium, and specific custimisations."""
    
    essential_rx = ['EX_ca2_e', 'EX_cl_e', 'EX_cobalt2_e', 'EX_cu2_e', 'EX_fe3_e', 'EX_fol_e', 
                    'EX_k_e', 'EX_mg2_e', 'EX_mn2_e', 'EX_pi_e', 'EX_so4_e', 'EX_ura_e', 'EX_zn2_e']
    
    full_medium = {rx:fill_value for rx in model.medium.keys()}

    for rx in essential_rx:
        full_medium[rx] = 10

    for rx in specific_reactions.keys():
        full_medium[rx] = specific_reactions[rx]
    
    return full_medium


def print_production_rates(sol):
    """Returns some metrics from a cobra solution in an easily readable format."""

    production_reactions = ["EX_but_e", "EX_ac_e", "EX_etoh_e", "EX_btoh_e", "EX_acetone_e"]
    print('Production fluxes:')
    print(sol.fluxes[production_reactions])


def plot_production_stats(solution_df):
    """Makes a simple plot to visualise the output of get_productions()."""

    plot_df = solution_df.copy()
    plot_df["reaction"] = solution_df.index
    plot_df.rename(columns={"pFBA sol": "pFBA", "minimum": "FVA min", "maximum": "FVA max"}, inplace=True)
    plot_df = plot_df.melt(id_vars="reaction", var_name="solution type", value_name="flux")

    sns.scatterplot(data=plot_df, x="reaction", y="flux", hue="solution type")


def plot_flux_envelopes(model, reactions: list = None, medium: dict = None, BM_func: str = "Growth"):
    """Plots production envelopes for each of the liste reactions in a single figure.

    Args:
        model: cobrapy GEM
        reactions (list, optional): Reactions to use for plot production envelope for. If None, ABE production reactions are used. Defaults to None.
        medium (dict, optional): Medium to use. If None the model medium is used. Defaults to None.
        BM_func (str, optional): Biomass function id. Defaults to "Growth".
    """

    if reactions is None:
        reactions = ["EX_but_e", "EX_ac_e", "EX_etoh_e", "EX_btoh_e", "EX_acetone_e"]

    colors = ['b', 'g', 'r', 'c', 'm', 'y']

    with model:
        
        if medium is not None:
            model.medium = medium

        # create subplots with 3 columns for each row
        num_reactions = len(reactions)

        num_rows = num_reactions // 3 if num_reactions % 3 == 0 else (num_reactions // 3) + 1
        _, axes = plt.subplots(num_rows, 3, figsize=(10, 3 * num_rows))

        for i, (rx, color) in enumerate(zip(reactions, colors)):
            # calculate the current row and column index in the grid
            row_index = i // 3
            col_index = i % 3
            
            # get the current subplot
            ax = axes[row_index, col_index] if num_reactions > 3 else axes[i]
            
            prod_env = production_envelope(model, [BM_func], objective=rx)

            # plot on the current subplot with the specified color
            prod_env.plot(kind='line', y=['flux_maximum', 'flux_minimum'], x=BM_func, ax=ax, color = color, legend=False)
            ax.fill_between(prod_env[BM_func], prod_env["flux_minimum"], prod_env["flux_maximum"], alpha=0.2, color=color)
            
            # set title for the subplot
            ax.set_title(rx)

        # remove empty subplots, if any
        for i in range(num_reactions, num_rows * 3):
            if num_reactions > 3:
                row, col = divmod(i, 3)
                axes[row, col].axis('off')
            else:
                axes[i].axis('off')

    plt.tight_layout()


def plot_flux_ranges(model = None, medium: dict=None, reactions: list=None, fva_sol: pd.DataFrame = None, 
                     log_scale: bool = True, loopless: bool = True, fraction_of_optimum: float = 1.0):
    """Plot flux ranges of reactions. Either takes in an existing FVA solution, or calculates a new one.

    Args:
        model (cobrapy model, optional): genome scale model. Defaults to None.
        medium (dict, optional): Medium. If None, default medium is used. Defaults to None.
        reactions (list, optional): reactions to plot flux ranges for. Defaults to None.
        fva_sol (pd.DataFrame, optional): Existing FVA solution. If given, new solution is not calculated. Defaults to None.
        log_scale (bool, optional): Symetric log-scale on x-axis. Defaults to True.
        loopless (bool, optional): use loopless FVA for analysis. Defaults to True.
    """

    if model is None and fva_sol is None:
        raise ValueError("Either model or existing FVA solution must be provided.")
    
    if model is not None and reactions is None:
        raise ValueError("Reactions must be provided if solution is to be calculated.")
    
    # calculate FVA result
    if fva_sol is None:
        with model:
            if medium is not None:
                model.medium = medium

            fva_sol = flux_variability_analysis(model, reactions, loopless=loopless, fraction_of_optimum=fraction_of_optimum)
    
    # plot figure
    
    plt.figure()
    
    for idx, row in fva_sol.iterrows():
        min_val = row['minimum']
        max_val = row['maximum']

        plt.plot([min_val, max_val], [idx, idx], marker='o', linestyle='-')

    if log_scale:
        plt.xscale('symlog')
        
    plt.xlabel('Flux (mmol/gDW/h)')
    plt.ylabel('Reaction')