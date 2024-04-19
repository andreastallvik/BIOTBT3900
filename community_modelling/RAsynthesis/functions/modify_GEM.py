"""Adapted from reframed to cobrapy"""

from cobra import Metabolite


def add_ratio_constraint_cobra(model, r_num, r_den, ratio):
        """ Add a flux ratio constraint to the model. Adapted from reframed to cobrapy by Andrea Stallvik.

        Arguments:
            r_num (str): id of the numerator
            r_den (str): id of the denominator
            ratio (float): ratio value

        Returns:
            str : identifier of the pseudo-metabolite
        """

        if r_num not in model.reactions or r_den not in model.reactions:
            raise KeyError(f"Invalid reactions in ratio {r_num}/{r_den}")

        pseudo_c_id = "pseudo"
        pseudo_m_id = f"ratio_{r_num}_{r_den}"

        pseudo_m = Metabolite(pseudo_m_id, compartment=pseudo_c_id)

        model.reactions.get_by_id(r_num).add_metabolites({pseudo_m: 1})
        model.reactions.get_by_id(r_den).add_metabolites({pseudo_m: -ratio})

        return pseudo_m