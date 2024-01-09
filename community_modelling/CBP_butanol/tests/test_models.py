"""Model validation tests."""

import pytest
import cobra


def validate_pathway(model, A, B, A_2=None) -> bool:
    """Tests whether the model network supports production of metabolite B from metabolite A."""
    
    test_model = model.copy()

    test_model.add_boundary(test_model.metabolites.get_by_id(A), type="sink", reaction_id="SK_A")
    if A_2 != None: # adds the option to add a second metabolite
        test_model.add_boundary(test_model.metabolites.get_by_id(A_2), type="sink", reaction_id="SK_A_2")

    test_model.add_boundary(test_model.metabolites.get_by_id(B), type="sink", reaction_id="DM_B")

    test_model.objective="DM_B"

    sol = test_model.slim_optimize()
    return sol > 0 


@pytest.fixture(params=["GEMs/NJ4_curated.xml"])#params=["GEMs/NJ4.xml", "GEMs/M5.xml"]
def model(request):
    return cobra.io.read_sbml_model(request.param)


def test_ATP_generation_cycle(model):
    """Tests whether model contains ATP-generation cycle."""
    
    test_model = model.copy()

    # set the lower bound of ATMP function to 0
    test_model.reactions.ATPM.lower_bound = 0

    # set the objective function to max ATP maintance
    test_model.objective = "ATPM"

    # set lower bound of all exhanges to 0 (no uptake allowed)
    for reaction in test_model.reactions:
            if reaction.boundary:
                    reaction.lower_bound = 0

    try:
        pfba_solution = cobra.flux_analysis.pfba(test_model)
        assert pfba_solution.fluxes["ATPM"] == 0, "Energy-generating cycle detected (ATP)."

    except cobra.exceptions.Infeasible:
        assert False, "Infeasible solution found."


def test_NADH_generation_cycle(model):
    """Tests whether model contains NADH-generation cycle."""
    
    test_model = model.copy()

    # set the lower bound of ATMP function to 0
    test_model.reactions.ATPM.lower_bound = 0

    # add a demand reaction for NADH
    test_model.add_boundary(test_model.metabolites.nadh_c, type="deman", reaction_id="DM_nadh_c")

    # set the objective function to be consumption of NADH
    test_model.objective = "DM_nadh_c"

    # set lower bound of all exhanges to 0 (no uptake allowed)
    for reaction in test_model.reactions:
            if reaction.boundary:
                    reaction.lower_bound = 0

    try:
        pfba_solution = cobra.flux_analysis.pfba(test_model)
        assert pfba_solution.fluxes["DM_nadh_c"] == 0, "Energy-generating cycle detected (NADH)."

    except cobra.exceptions.Infeasible:
        assert False, "Infeasible solution found."


def test_pathways(model):
    """Test some key conversions in the ABE fermentation pathway."""

    assert validate_pathway(model, "xyl__D_c", "accoa_c"), "Pathway xyl__D_c -> accoa_c not feasible."
    assert validate_pathway(model, "accoa_c", "etoh_c"), "Pathway accoa_c -> etoh_c not feasible."
    assert validate_pathway(model, "accoa_c", "ac_c"), "Pathway accoa_c -> ac_c not feasible."
    assert validate_pathway(model, "accoa_c", "acetone_c"), "Pathway accoa_c -> acetone_c not feasible."
    assert validate_pathway(model, "accoa_c", "but_c"), "Pathway accoa_c -> but_c not feasible."
    assert validate_pathway(model, "accoa_c", "btoh_c"), "Pathway accoa_c -> btoh_c not feasible."
    assert validate_pathway(model, A= "ac_c", A_2 = "but_c", B="btoh_c"), "Pathway ac_c + but_c -> btoh_c not feasible."