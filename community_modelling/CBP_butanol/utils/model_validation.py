"""
Helper functions for checking validity of automatically reconstructed models.
"""

import cobra
from cobra.io import read_sbml_model
from cobra.exceptions import Infeasible
import traceback
import networkx as nx
import pyvis.network


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


def get_exhange(reactions):
    """Returns first exhange reaction in list of reactions, returns false if no exchange reaction is found."""
    for reaction in reactions:
        if reaction.id.startswith("EX_"):
            return reaction
    return False


def add_exhanges_model(test_model, A):
    """Adds exhange reactions for a metabolite A to the model. We assume that the metabolite is in the cytosol."""
    
    if test_model.metabolites.get_by_id(A).compartment != "C_e":
        # create metabolite A_e
        A_e = cobra.Metabolite(A + "_e", compartment = "C_e")
        # add a reaction converting A to A_e
        transport_A = cobra.Reaction('transport_' + A)
        transport_A.name = 'transport_A'
        transport_A.add_metabolites({
            test_model.metabolites.get_by_id(A): -1.0,
            A_e: 1.0,
        })
        
        test_model.add_reactions([transport_A])

        # create an exhange reaction for A_e

        test_model.add_boundary(test_model.metabolites.get_by_id(A_e.id), type="exchange")
        
        # get the name of the exhange reaction
        adjacent_rx = [rx for rx in test_model.metabolites.get_by_id(A_e.id).reactions]
        A_exchange = get_exhange(adjacent_rx)
    
    else:
        print("Error: metabolite is extracellular!")
        # A_e = test_model.metabolites.get_by_id(A)

        # # get the name of the exhange reaction
        # adjacent_rx = [rx for rx in test_model.metabolites.get_by_id(A_e.id).reactions]
        # A_exchange = get_exhange(adjacent_rx)

        # if A_exchange == False:
        #     print("Likely an error in the model, the should be exhange rx for external metabolite.")
        
    return test_model, A_exchange


def exists_path(model, A, B) -> bool:
    """Tests that there exists a path from metabolite A to metabolite B in the model."""
    
    test_model = model.copy()

    if A not in test_model.metabolites:
        print("Metabolite " + A + " not in model.")
        return False
    elif B not in test_model.metabolites:
        print("Metabolite " + B + " not in model.")
        return False
    
    test_model, A_exchange = add_exhanges_model(test_model, A)
    test_model, B_exchange = add_exhanges_model(test_model, B)

    # set uptake of the exhange to be -10
    test_model.reactions.get_by_id(A_exchange.id).lower_bound = -10
    
    # set production of B as objective
    test_model.objective = B_exchange.id

    # set ATPM to 0
    test_model.reactions.get_by_id("ATPM").lower_bound = 0

    sol = test_model.optimize()
    if sol.objective_value > 0:
        return sol.fluxes[sol.fluxes > 0]
    else:
        return False


def production_possible(model, reaction, medium = None) -> bool:
    """Tests whether the model is able to carry flux through the reaction at steady-state on the specified medium."""
    
    test_model = model.copy()

    if reaction not in test_model.reactions:
        print("Reaction " + reaction + " not in model.")
        return False

    if medium is None:
        medium = test_model.medium 

    test_model.medium = medium

    test_model.objective = reaction
    
    sol = test_model.optimize()
    
    return sol.objective_value > 0


def reactions_exist(model, reactions) -> bool:
    """Tests that the model contains the specified reactions."""
    # TODO: implement
    pass


def create_graph(model, metabolites = None):
    """Create a (bipartite) graph representation of the metabolic network by iterating over each reaction."""

    G = nx.Graph()

    if metabolites is not None:
        set_of_reactions = set()
        for m in metabolites:
            set_of_reactions.add(model.metabolites.get_by_id(m).reactions)

        reactions = {element for frozenset in set_of_reactions for element in frozenset}
    else:
        reactions = model.reactions

    for reaction in reactions:

        reaction_dict = reaction.metabolites

        consumed = [key.id for key, value in reaction_dict.items() if value < 0]
        produced = [key.id for key, value in reaction_dict.items() if value > 0]

        G.add_nodes_from(consumed, bipartite=0)
        G.add_nodes_from(produced, bipartite=0)
        G.add_node(reaction.id, bipartite=1, color = "red", shape = "square")

        for i in consumed:
            G.add_edge(i, reaction.id)

        for j in produced:
            G.add_edge(reaction.id, j)

    return G


def visualise_graph(G, filename):
    """Visualise a graph using pyvis, write as html file."""
    nt = pyvis.network.Network('1000px', '1800px')
    nt.from_nx(G)
    filepath = filename + ".html"
    nt.show(filepath, notebook=False)


# main - to troubleshoot

# nj4 = read_sbml_model('community_modelling/CBP_butanol/GEMs/NJ4.xml')
# print(exists_path(nj4, "glc__D_c", "ac_c"))