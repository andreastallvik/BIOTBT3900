"""
Helper functions for checking validity of automatically reconstructed models.
"""

import cobra
from cobra.io import read_sbml_model
from cobra.exceptions import Infeasible
import traceback


def energy_generation_cycle(model) -> bool:
    """Tests for energy-generation cycles in the model."""
    
    test_model = model.copy()

    # set the lower bound of ATMP function to 0
    test_model.reactions.ATPM.lower_bound = 0

    # set the objective function to max ATP maintance
    test_model.objective = "ATPM"

    # set lower bound of all exhanges to 0 (no uptake allowed)
    for reaction in test_model.reactions:
            if reaction.boundary:
                    reaction.lower_bound = 0

    # solve for max ATPM reaction with pFBA
    try:
        pfba_solution = cobra.flux_analysis.pfba(test_model)
        print("ATPM flux: " + str(pfba_solution.fluxes["ATPM"]))
        if pfba_solution.fluxes["ATPM"] > 0:
            print("Energy-generating cycle detected.")
            return True
        else:
            return False
    except Infeasible:
        traceback.print_exc()


def growth_possible(model, medium = None) -> bool:
    """Tests that the model can grow on the specified medium."""

    test_model = model.copy()

    if medium is None:
        test_model.medium = medium
    
    sol = test_model.optimize()
    
    if sol.objective_value == 0:
        return False
    else:
        return True


def production_possible(model, metabolites, medium = None) -> bool:
    """Tests that the model can produce the specified metabolite on the specified medium."""
    
    test_model = model.copy()

    if medium is None:
        test_model.medium = medium

    sol = test_model.optimize()



def reactions_exist(model, reactions) -> bool:
    """Tests that the model contains the specified reactions."""
    # TODO: implement
    pass