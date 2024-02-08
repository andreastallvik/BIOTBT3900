"""Preform all modifications to the out-of-the-box carveme reconstruction of NJ4 and save the modified model as NJ4_curated"""

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from config import ROOT_DIR

print("Modifying NJ4 ...")

# load nj4
print("loading nj4 model...")
nj4 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "NJ4.xml"))

# load universal model
print("loading universal model...")
universal_model = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

print("retrieving reactions...")
# butabol transport reaction
BTOHt = universal_model.reactions.get_by_id('BTOHt')

# butyrate transport reaction
BUTt = universal_model.reactions.get_by_id('BUTt')

# acetoacetate -> acetone + CO2
ADCi = universal_model.reactions.get_by_id("ADCi")

# acetone transport reaction
ACEt = universal_model.reactions.get_by_id("ACEt")

# acetate transport reaction
ACtr = universal_model.reactions.get_by_id("ACtr")

# Acyl-CoA dehydrogenase (butanoyl-CoA) - added due to an article claiming this reaction is solely NADH dependant. consider incorporating ferrodoxin
ACOAD1 = universal_model.reactions.get_by_id("ACOAD1")

# ferrodoxin oxidoreductase reaction
POR_syn = universal_model.reactions.get_by_id("POR_syn")

# Ferredoxin:NAD+ reductase
FNRR = universal_model.reactions.get_by_id("FNRR")

# (FeFe)-hydrogenase, EC 1.12.1.4, KeGG rx: R09508
HYDA = cobra.Reaction('HYDA')
HYDA.add_metabolites({
    nj4.metabolites.get_by_id("h2_c"): -2.0,
    nj4.metabolites.get_by_id("nad_c"): -1.0,
    nj4.metabolites.get_by_id("fdxo_2_2_c"): -2.0,
    nj4.metabolites.get_by_id("h_c"): 5.0,
    nj4.metabolites.get_by_id("nadh_c"): 1.0,
    nj4.metabolites.get_by_id("fdxrd_c"): 2.0,
})
HYDA.bounds = (-1000, 1000)

# make these reaction reversible like god intended
# CoA- transferase reactions by enzyme 2.8.3.8
nj4.reactions.get_by_id("BUTCT").bounds = (-1000, 1000) 
nj4.reactions.get_by_id("ACACCT").bounds = (-1000, 1000) 
nj4.reactions.get_by_id("BUTCT2").bounds = (-1000, 1000) 


print("adding reactions...")
nj4.add_reactions([BTOHt, BUTt, ADCi, ACEt, ACtr, ACOAD1, POR_syn, FNRR, HYDA])
nj4.add_boundary(nj4.metabolites.get_by_id('btoh_e'), type='exchange', reaction_id='EX_btoh_e', lb=0, ub=1000)
nj4.add_boundary(nj4.metabolites.get_by_id('acetone_e'), type='exchange', reaction_id='EX_acetone_e', lb=0, ub=1000)

print("writing model as NJ4_curated...")
write_sbml_model(nj4, str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "NJ4_curated.xml"))