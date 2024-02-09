import cobra
from cobra.io import read_sbml_model, write_sbml_model
from config import ROOT_DIR

print("Modifying M5 ...")

# load nj4
print("loading m5 model...")
m5 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "M5.xml"))

# load universal model
print("loading universal model...")
universal_model = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

print("retrieving reactions...")

### -------- xylan degradation reactions -------- ###

# xyl4_e metabolite
xyl4_e = cobra.Metabolite(
    'xyl4_e',
    formula='C20H34O17',
    name='Xylotetraose',
    compartment='C_e')

# glcur_e <-> glcur_c transport rx
glcur_transport = cobra.Reaction('glcur_transport')
glcur_transport.add_metabolites({
    m5.metabolites.get_by_id("glcur_c"): -1.0,
    m5.metabolites.get_by_id("glcur_e"): 1.0,
})
glcur_transport.bounds = (-1000, 1000)

# GLCURS1 rx, extracellular edition
GLCURS1_e = cobra.Reaction('GLCURS1_e')
GLCURS1_e.add_metabolites({
    m5.metabolites.get_by_id("h2o_e"): -1.0,
    m5.metabolites.get_by_id("xylan4_e"): -1.0,
    m5.metabolites.get_by_id("glcur_e"): 1.0,
    xyl4_e: 1.0,
})
GLCURS1_e.bounds = (-1000, 1000)

# XYLOS1 rx, extracellular edition
XYLOS1_e = cobra.Reaction('XYLOS1_e')
XYLOS1_e.add_metabolites({
    m5.metabolites.get_by_id("h2o_e"): -3.0,
    xyl4_e: -1.0,
    m5.metabolites.get_by_id("xyl__D_e"): 4.0,
})
XYLOS1_e.bounds = (-1000, 1000)

### -------- ABE fermentation reactions -------- ###
# butabol transport reaction
BTOHt = universal_model.reactions.get_by_id('BTOHt')

# butyrate transport reaction
BUTt = universal_model.reactions.get_by_id('BUTt')

# acetate transport reaction
ACtr = universal_model.reactions.get_by_id("ACtr")

# ethanol transport reaction
ETOHt = universal_model.reactions.get_by_id("ETOHt")

# Ferredoxin:NAD+ reductase
FNRR = universal_model.reactions.get_by_id("FNRR")

# Acyl-CoA dehydrogenase (butanoyl-CoA) - added due to an article claiming this reaction is solely NADH dependant. consider incorporating ferrodoxin
ACOAD1 = universal_model.reactions.get_by_id("ACOAD1")

# block reaction as there is no evidence for existance of this enzyme (2.8.3.8) in the bacteria
m5.reactions.get_by_id("BUTCT").bounds = (0, 0)

# (FeFe)-hydrogenase, EC 1.12.1.4, KeGG rx: R09508 - evidence for this in prev paper from same authors
HYDA = cobra.Reaction('HYDA')
HYDA.add_metabolites({
    m5.metabolites.get_by_id("h2_c"): -2.0,
    m5.metabolites.get_by_id("nad_c"): -1.0,
    m5.metabolites.get_by_id("fdxo_2_2_c"): -2.0,
    m5.metabolites.get_by_id("h_c"): 5.0,
    m5.metabolites.get_by_id("nadh_c"): 1.0,
    m5.metabolites.get_by_id("fdxrd_c"): 2.0,
})
HYDA.bounds = (-1000, 1000)

print("adding reactions...")
m5.add_reactions([glcur_transport, GLCURS1_e, XYLOS1_e, BTOHt, BUTt, ACtr, ETOHt, FNRR, ACOAD1, HYDA])
m5.add_boundary(m5.metabolites.get_by_id('btoh_e'), type='exchange', reaction_id='EX_btoh_e', lb=0, ub=1000)

print("writing model as M5_curated...")
write_sbml_model(m5, str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "M5_curated.xml"))