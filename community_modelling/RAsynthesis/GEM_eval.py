"""
Functions to evaluate GEMs.
"""

from cobra import flux_analysis
from cobra.exceptions import Infeasible
from warnings import warn

def energy_generating_cycle_test(model):
    """Testing for energy-generating cycles"""

    test_model = model.copy()

    #set lower bound of ATMP function to 0
    test_model.reactions.ATPM.lower_bound = 0

    #Set the objective function to max ATP maintance
    test_model.objective = "ATPM"

    #Set lower bound of all exhanges to 0
    for reaction in test_model.reactions:
            if reaction.boundary:
                    reaction.lower_bound = 0

    #Maximise the ATPM reaction using pFBA
    try:
        pfba_solution = flux_analysis.pfba(test_model)
        print("ATPM flux: " + str(pfba_solution.fluxes["ATPM"]))
        return pfba_solution
    except Infeasible:
        print("Infesible solution")

def saa_loop(sol, community = True):
    """Tests that the split patway in SAA module is not behaving like a loop"""

    if community:
        reactions = ["R_HPPHD_SAL9", "R_DHPPSA_SAL9", "R_DLDH_SAL9", "R_HPLSA_SAL9"]

    else:
        reactions = ["R_HPPHD", "R_DHPPSA", "R_DLDH", "R_HPLSA"]

    is_loop = False

    for reaction in reactions:
        if (reaction == "R_DLDH_SAL9") or (reaction == "R_DLDH"):
            if sol.values[reaction] > 0:
                warn(reaction + " has flux value of " +str(sol.values[reaction]))
                is_loop = True
        else:
            if sol.values[reaction] < 0:
                warn(reaction + " has flux value of " +str(sol.values[reaction]))
                is_loop = True
    
    return is_loop