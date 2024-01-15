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

# xyl_e <-> xyl_c transport rx
xyl4_transport = cobra.Reaction('xyl4_transport')
xyl4_transport.add_metabolites({
    m5.metabolites.get_by_id("xyl4_c"): -1.0,
    xyl4_e: 1.0
})

# glcur_e <-> glcur_c transport rx
glcur_transport = cobra.Reaction('glcur_transport')
glcur_transport.add_metabolites({
    m5.metabolites.get_by_id("glcur_c"): -1.0,
    m5.metabolites.get_by_id("glcur_e"): -1.0,
})

# GLCURS1 rx, extracellular edition
GLCURS1_e = cobra.Reaction('GLCURS1_e')
GLCURS1_e.add_metabolites({
    m5.metabolites.get_by_id("h2o_e"): -1.0,
    m5.metabolites.get_by_id("xylan4_e"): -1.0,
    m5.metabolites.get_by_id("glcur_e"): 1.0,
    xyl4_e: 1.0,
})

# XYLOS1 rx, extracellular edition
XYLOS1_e = cobra.Reaction('XYLOS1_e')
XYLOS1_e.add_metabolites({
    m5.metabolites.get_by_id("h2o_e"): -3.0,
    xyl4_e: -1.0,
    m5.metabolites.get_by_id("xyl__D_e"): 4.0,
})

### -------- ABE fermentation reactions -------- ###
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

# ethanol transport reaction
ETOHt = universal_model.reactions.get_by_id("ETOHt")

# CoA- transferase reactions by enzyme 2.8.3.8
ACACCT = universal_model.reactions.get_by_id("ACACCT")
BUTCT2 = universal_model.reactions.get_by_id("BUTCT2")

# make these reaction reversible like god intended
m5.reactions.get_by_id("BUTCT").bounds = (-1000, 1000) 

print("adding reactions...")
m5.add_reactions([xyl4_transport, glcur_transport, GLCURS1_e, XYLOS1_e, BTOHt, BUTt, ADCi, ACEt, ACtr, ETOHt, ACACCT, BUTCT2])
m5.add_boundary(m5.metabolites.get_by_id('btoh_e'), type='exchange', reaction_id='EX_btoh_e', lb=0, ub=1000)
m5.add_boundary(m5.metabolites.get_by_id('acetone_e'), type='exchange', reaction_id='EX_acetone_e', lb=0, ub=1000)

print("writing model as M5_curated...")
write_sbml_model(m5, str(ROOT_DIR / "community_modelling" / "CBP_butanol" / "GEMs" / "M5_curated.xml"))