"""
Helper functions for preforming/displaying static simulations of the CBP butanol model.
"""

from cobra.flux_analysis import pfba, flux_variability_analysis
import pandas as pd


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
    
    
def get_specific_medium(model, specific_reactions: list):
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