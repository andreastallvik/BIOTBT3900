"""Preform all modifications to the out-of-the-box carveme reconstruction of NJ4 and save the modified model as NJ4_curated"""

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

# acetoacetyl-CoA + butyrate -> acetoacetate + butanoyl-CoA
nj4.reactions.get_by_id("BUTCT2").bounds = (-1000, 1000) # make this reaction reversible like god intended

print("adding reactions...")
nj4.add_reactions([BTOHt, BUTt, ADCi, ACEt, ACtr])
nj4.add_boundary(nj4.metabolites.get_by_id('btoh_e'), type='exchange', reaction_id='EX_btoh_e')
nj4.add_boundary(nj4.metabolites.get_by_id('acetone_e'), type='exchange', reaction_id='EX_acetone_e')

print("writing model as NJ4_curated...")
write_sbml_model(nj4, str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "NJ4_curated.xml"))