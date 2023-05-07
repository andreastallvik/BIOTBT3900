from cobra import flux_analysis
from cobra.exceptions import Infeasible

"""Testing for energy-generating cycles"""
def energy_generating_cycle_test(model):

    test_model = model.copy()

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