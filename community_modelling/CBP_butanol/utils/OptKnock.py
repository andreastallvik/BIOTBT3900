import cobra
import numpy as np
import pandas as pd
from itertools import combinations
import copy


def phonyOptKnock(model, target, reactions=None, knockouts=1):
    """ Phony implementation of OptKnock to identify mutants 
    with growth-coupled target metabolite production.

    Parameters
    ----------
    model : cobra model object
            Model object.
    
    target : string
             Reaction id of target reaction to be improved through model knockouts.

    reactions : cobra reaction object(s) 
                Reactions to be considered for knockout. By default
                all model reactions.

    knockouts : int
                Number of knockouts to consider, single (1),
                double (2), triple (3), etc.

    Returns
    ----------
    results : pandas DataFrame object
              Results for all proposed mutants containing the columns:

              reactions: knocked-out reaction(s)

              target: flux of target reaction when optimizing the model objective function

              biomass: corresponding optimal objective value of the model (i.e. growth if biomass is the objective)

              fva_min: minimal flux through the target reaction at optimal growth

              fva_max: maximal flux through the target reaction at optimal growth
    """

    del_model = copy.deepcopy(model)

    # By default all model reactions
    if reactions == None:
        reactions = del_model.reactions

    # Remove biomass, exchange, and reactions without associated genes
    reactions.remove(list(cobra.util.solver.linear_reaction_coefficients(del_model).keys())[0])
    reactions = [rxn for rxn in reactions if rxn not in cobra.medium.find_boundary_types(del_model, "exchange")]
    reactions = [rxn for rxn in reactions if not "s0001" in rxn.gene_reaction_rule]

    # All unique knockout combinations
    if knockouts > 1:
        candidates = [i for i in combinations(reactions, knockouts)]
    else:
        candidates = list(reactions)
    
    result = pd.DataFrame(index=range(len(candidates)),columns=["reactions", "target", "biomass", "fva_min", "fva_max"])

    sol_wt = del_model.optimize()

    # Constrain and simulate fluxes
    for i, cand in enumerate(candidates):
        skip = False
        if knockouts == 1:
            cand = [cand]
        else:
            cand = list(cand)

        old_lower_bounds = np.zeros(knockouts)
        old_upper_bounds = np.zeros(knockouts)
        for j, c in enumerate(cand):
            old_lower_bounds[j] = cand[j].lower_bound
            old_upper_bounds[j] = cand[j].upper_bound
            c.lower_bound = 0.0
            c.upper_bound = 0.0

        # Simulate fluxes
        sol = del_model.optimize()

        # Skip if growth is below 10% of wild type
        if sol.objective_value < 0.1 * sol_wt.objective_value:
            skip = True

        # Simulate flux variability, skip if infeasible
        try:
            sol_fva = cobra.flux_analysis.flux_variability_analysis(del_model, target)
        except:
            skip = True

        if not skip:
            # Add results to dataframe
            result["reactions"][i] = [c.id for c in cand]
            result["target"][i] = sol.fluxes[target]
            result["biomass"][i] = sol.objective_value
            result["fva_min"][i] = sol_fva["minimum"][0]
            if sol_fva["maximum"][0] < 1e-10:
                result["fva_max"][i] = 0
            else:
                result["fva_max"][i] = sol_fva["maximum"][0]
                

        # Reset previous flux bounds
        for j, rxn in enumerate(cand):
            rxn.upper_bound = old_upper_bounds[j]
            rxn.lower_bound = old_lower_bounds[j]
    
    # Post-process results - remove nan rows and sort
    result = result.dropna()
    result = result.sort_values(by=['fva_min'], ascending=[False])
    return result
