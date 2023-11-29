"""
Script for running 16 experiments with all inoculation ratio combinations for CAL2:SAL9:MAM2,
in order to re-create parts of fig. 4 D from the RA paper.
Saves result as csv file in /results folder.
NOTE: Filepaths are hardocded, so scripts needs to be run from within the /function directory for anything to work.
"""

import itertools
import pandas as pd
from tqdm import tqdm
import data_analysis
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_glc_triculture


## load and prepare models
print("loading exp. data")

# load exp. data yields, but magnified
relative_abundance_df = pd.read_csv("../exp_data/subpop_data_glc.csv")
products_df = pd.read_csv("../exp_data/conc_data_glc.csv")
od600_df = pd.read_csv("../exp_data/od600_glc.csv")
subpop_df, conc_df, OD_df = data_analysis.process_data(relative_abundance_df, products_df, od600_df)
CA_yield, SAA_yield, RA_yield = data_analysis.get_yields_glc_xyl(conc_df)

print("loading and preparing GEMs")

# ----- read cobrapy models -----

CAL2_cobra = read_sbml_model("../GEMs/CAL2.xml")
SAL9_cobra = read_sbml_model("../GEMs/SAL9.xml")
MAM2_cobra = read_sbml_model("../GEMs/MAM2.xml")

# ----- add yield ratio constraint -----

# knock out GLCtex_copy1 in SAL9, just so I have one less thing to worry about...
CAL2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
SAL9_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
MAM2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

add_ratio_constraint_cobra(CAL2_cobra, "34DHCINMt", "GLCtex_copy2", CA_yield);
add_ratio_constraint_cobra(SAL9_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
add_ratio_constraint_cobra(MAM2_cobra, "RAt", "GLCtex_copy2", RA_yield);

inocculation_ratios = [
    (1,1,2),
    (1,2,1),
    (2,1,1),
    (1,1,1),
    (1,2,2),
    (2,1,2),
    (2,2,1),
    (1,2,3),
    (1,3,2),
    (2,1,3),
    (2,3,1),
    (3,1,2),
    (3,2,1),
    (1,3,3),
    (3,1,3),
    (3,3,1),
]


print("setting up experiment")


tot_BM_list = []
tot_RAs = []
inoc_ratio_list = []
glc_xyl_ratio_list = []

for inocculation_ratio in tqdm(inocculation_ratios):
    
    try:

        sim = simulate_glc_triculture(CAL2_cobra, SAL9_cobra, MAM2_cobra, 
                                          initial_pop=2.7e-3, initial_pop_ratio=inocculation_ratio, 
                                          adjust_atp_requirements=True, km_glc_adj=1000, glc_vmax_adjustment=0.45)

        tot_BM = sum(sim.total_biomass.drop(columns=["cycle"], inplace=False).iloc[-1])
        tot_RA = sim.get_metabolite_time_series()["rosma_e"].iloc[-1]

        tot_BM_list.append(tot_BM)
        tot_RAs.append(tot_RA)
        inoc_ratio_list.append(inocculation_ratio)

    except Exception as e:
        print("not able to process:", inocculation_ratio)
        print(f"An exception occurred: {e}")

results_df = pd.DataFrame({'inoculation_ratio': inoc_ratio_list, 'total_biomass': tot_BM_list, 'total_RA': tot_RAs})

results_df.to_csv("../results/fig_4D_29nov", index=False)
