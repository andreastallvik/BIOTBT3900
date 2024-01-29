from cobra.io import read_sbml_model
from utils.OptKnock import phonyOptKnock
from utils import static_sim

nj4 = read_sbml_model("GEMs/NJ4_curated.xml")

medium = static_sim.get_specific_medium(nj4, {"EX_xyl__D_e": 10})

with nj4:

    nj4.medium = medium

    optKnock = phonyOptKnock(nj4, "EX_but_e", knockouts=1)

    optKnockdf = optKnock[optKnock['fva_min'] > 0]

    print(len(optKnockdf))

    print(optKnockdf.sort_values('fva_min', ascending=False))