from cobra.io import read_sbml_model
from utils.OptKnock import phonyOptKnock
from utils import static_sim

m5 = read_sbml_model("GEMs/M5_curated.xml")

medium = static_sim.get_specific_medium(m5, {"EX_xylan4_e": 1, "EX_xylan8_e": 0.5, "EX_nh4_e": 10})

# added constraint to set a direction
m5.reactions.HACD1i.bounds = (-1000, 0)
m5.reactions.HBCO_nadp.bounds = (0, 1000)
m5.reactions.ALCD2x.bounds = (-1000, 0)
m5.reactions.ALCD2y.bounds = (-1000, 0)

with m5:

    m5.medium = medium

    optKnock = phonyOptKnock(m5, "EX_but_e", knockouts=1)

    optKnockdf = optKnock[optKnock['fva_min'] > 0]

    print(len(optKnockdf))

    print(optKnockdf.sort_values('fva_min', ascending=False))