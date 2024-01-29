from cobra.io import read_sbml_model
from utils.OptKnock import phonyOptKnock

nj4 = read_sbml_model("GEMs/NJ4_curated.xml")

optKnock = phonyOptKnock(nj4, "EX_ac_e", knockouts=1)

optKnockdf = optKnock[optKnock['fva_min'] > 0]

print(len(optKnockdf))

print(optKnockdf.sort_values('fva_min', ascending=False))