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
    
    
def get_specific_medium(model, specific_reactions: dict):
    """Get a dictionary of max uptake rates for a model, with 0.1 for each reaction in the original medium, and specific custimisations."""
    
    essential_rx = ['EX_ca2_e', 'EX_cl_e', 'EX_cobalt2_e', 'EX_cu2_e', 'EX_fe3_e', 'EX_fol_e', 
                    'EX_k_e', 'EX_mg2_e', 'EX_mn2_e', 'EX_pi_e', 'EX_so4_e', 'EX_ura_e', 'EX_zn2_e']
    
    full_medium = {rx:0.1 for rx in model.medium.keys()}

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


def plot_flux_envelopes(model, reactions: list = None, medium: dict = None):

    if reactions is None:
        reactions = ["EX_but_e", "EX_ac_e", "EX_etoh_e", "EX_btoh_e", "EX_acetone_e"]

    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    
    with model:
        
        if medium is not None:
            model.medium = medium

        # create subplots with 3 columns for each row
        num_reactions = len(reactions)
        num_rows = num_reactions // 3 if num_reactions % 3 == 0 else (num_reactions // 3) + 1
        fig, axes = plt.subplots(num_rows, 3, figsize=(10, 3 * num_rows))

        for i, (rx, color) in enumerate(zip(reactions, colors)):
            # calculate the current row and column index in the grid
            row_index = i // 3
            col_index = i % 3
            
            # get the current subplot
            ax = axes[row_index, col_index] if num_reactions > 1 else axes
            
            prod_env = production_envelope(model, [rx])

            # plot on the current subplot with the specified color
            prod_env.plot(kind='line', x='flux_maximum', y=rx, xlabel="Growth rate", ax=ax, color=color)
            ax.fill_between(prod_env["flux_maximum"], prod_env[rx], alpha=0.2, color=color)
            
            # set title for the subplot
            ax.set_title(rx)

    plt.tight_layout()
