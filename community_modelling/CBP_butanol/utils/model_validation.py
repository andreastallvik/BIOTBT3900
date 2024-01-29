"""
Helper functions for checking validity of automatically reconstructed models.
"""

import cobra
from cobra.io import read_sbml_model
from cobra.exceptions import Infeasible
from cobra.flux_analysis import gapfill
import traceback
import networkx as nx
import pyvis.network


def ATP_generation_cycle(model):
    """Tests whether model contains ATP-generation cycle."""
    
    with model:

        # set the lower bound of ATMP function to 0
        model.reactions.ATPM.lower_bound = 0

        # set the objective function to max ATP maintance
        model.objective = "ATPM"

        # set lower bound of all exhanges to 0 (no uptake allowed)
        for reaction in model.reactions:
                if reaction.boundary:
                        reaction.lower_bound = 0

        # solve for max ATPM reaction with pFBA
        try:
            pfba_solution = cobra.flux_analysis.pfba(model)
            #print("ATPM flux: " + str(pfba_solution.fluxes["ATPM"]))
            if pfba_solution.fluxes["ATPM"] > 0:
                print("Energy-generating cycle detected.")
                return True
            else:
                return False
        except Infeasible:
            traceback.print_exc()


def NADH_generation_cycle(model):
    """Tests whether model contains NADH-generation cycle."""
    
    with model:

        # set the lower bound of ATMP function to 0
        model.reactions.ATPM.lower_bound = 0

        # add a demand reaction for NADH
        model.add_boundary(model.metabolites.nadh_c, type="deman", reaction_id="DM_nadh_c")

        # set the objective function to be consumption of NADH
        model.objective = "DM_nadh_c"

        # set lower bound of all exhanges to 0 (no uptake allowed)
        for reaction in model.reactions:
                if reaction.boundary:
                        reaction.lower_bound = 0

        try:
            pfba_solution = cobra.flux_analysis.pfba(model)
            if pfba_solution.fluxes["DM_nadh_c"] > 0:
                print("Energy-generating cycle detected.")
                return True
            else:
                return False
        except Infeasible:
            traceback.print_exc()
            

def energy_generation_cycle(model):
    """Tests whether model contains energy-generating cycle (ATP or NADH - returns false if neither is detected). """
    atp = ATP_generation_cycle(model)
    nadh = NADH_generation_cycle(model)
    return atp or nadh


def growth_possible(model, medium = None) -> bool:
    """Tests that the model can grow on the specified medium."""

    with model:

        if medium is not None:
            model.medium = medium
        
        sol = model.slim_optimize()
        
    if sol > 0:
        return True
    else:
        return False
    

def reactions_exist(model, reactions) -> bool:
    """Tests that the model contains all the reactions in the list."""
    
    missing_reactions = []

    with model:

        for reaction in reactions:
            if reaction not in model.reactions:
                missing_reactions.append(reaction)
    
    if missing_reactions:
        print("The following reactions are missing: " + str(missing_reactions))
        return False
    else:
        return True


def check_production(model, reactions, medium = None, fraction_of_optimum=1.0) -> bool:
    """Get the possible flux values for the reactions at any optimal solution."""
    
    with model:

        if not reactions_exist(model, reactions):
            raise ValueError("Input only reactions that are present in the model.")

        if medium is not None:
            model.medium = medium

        # run FVA for the reactions and check that both their scores are not 0
        fva_sol = cobra.flux_analysis.flux_variability_analysis(model, reactions, fraction_of_optimum=fraction_of_optimum)
    
    return fva_sol


def check_blocked_reactions(model, reactions, medium = None):
    """Checks whether any of the reactions are blocked reactions. 
    Returns the input reactions which are blocked. Uses cobra.flux_analysis.find_blocked_reactions, but allows for str id input on list of reactions."""
    
    with model:

        if medium is not None:
            model.medium = medium

        reaction_list = [model.reactions.get_by_id(r) for r in reactions]

        blocked_reactions = cobra.flux_analysis.find_blocked_reactions(model, reaction_list=reaction_list)

    if blocked_reactions:
        return blocked_reactions
    else:
        return False


def validate_pathway(model, A, B) -> bool:
    """Tests whether the model network supports production of metabolite B from metabolite A."""
    
    with model:    
        
        model.add_boundary(model.metabolites.get_by_id(A), type="sink", reaction_id="SK_A")
        model.add_boundary(model.metabolites.get_by_id(B), type="sink", reaction_id="DM_B")

        model.objective="DM_B"

        sol = model.slim_optimize()

    return sol > 0 


def find_reaction_connections(model, A, B):
    """Find the reaction(s) that connect metabolite A and B in a given model."""

    rx_A = model.metabolites.get_by_id(A).reactions
    rx_B = model.metabolites.get_by_id(B).reactions

    # get list of elements in both r1 and r2
    common = list(set(rx_A).intersection(rx_B))
    return [r.id for r in common]


def find_gapfilling_reactions(model, universal_model, reaction):
    """Using cobrapy gapfill function, find reactions from the universal model that when added will unblock reaction."""
    
    with model:
        model.objective = model.reactions.get_by_id(reaction)
        gapfilling_reactions = gapfill(model, universal_model)
        
    return [r.id for r in gapfilling_reactions[0]]


## archive(?) - need to figure out how to actually use these functions

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