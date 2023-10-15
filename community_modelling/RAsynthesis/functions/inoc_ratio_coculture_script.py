"""
Script for running 7 experiments with all inoculation ratio combinations for RAU2:RAD4,
in order to re-create parts of fig. 3 D from the RA paper.
Saves result as csv file in /results folder.
NOTE: Filepaths are hardocded, so scripts needs to be run from within the /function directory for anything to work.
"""

import itertools
import pandas as pd
from tqdm import tqdm
import data_analysis
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_coculture


## load and prepare models
print("loading exp. data")

# load exp. data yields, but magnified
CA_yield, SAA_yield, RA_yield = data_analysis.get_yields_coculture()
CA_yield = CA_yield * 10
SAA_yield = SAA_yield * 10
RA_yield = RA_yield * 10

print("loading and preparing GEMs")

# ----- read cobrapy models -----

RAU2_cobra = read_sbml_model("../GEMs/RAU2.xml")
RAD4_cobra = read_sbml_model("../GEMs/RAD4.xml")

# ----- add yield ratio constraint -----

# knock out GLCtex_copy1 in both models, just so I have one less thing to worry about...
RAU2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
RAD4_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

add_ratio_constraint_cobra(RAU2_cobra, "34DHCINMt", "GLCtex_copy2", CA_yield);
add_ratio_constraint_cobra(RAD4_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
add_ratio_constraint_cobra(RAD4_cobra, "RAt", "GLCtex_copy2", RA_yield);

inocculation_ratios = [
    (9, 1),
    (5, 1),
    (3, 1),
    (1, 1),
    (1, 3),
    (1, 5),
    (1, 9)
]

print("setting up experiment")


tot_BM_list = []
tot_CAs = []
tot_SAAs = []
tot_RAs = []
inoc_ratio_list = []

for inocculation_ratio in tqdm(inocculation_ratios):
    
    try:

        sim = simulate_coculture(RAU2_cobra, RAD4_cobra, initial_pop=2.e-3, 
                                 initial_pop_ratio=inocculation_ratio)

        tot_BM = sum(sim.total_biomass.drop(columns=["cycle"], inplace=False).iloc[-1])
        tot_CA = sim.get_metabolite_time_series()["34dhcinm_e"].iloc[-1]
        tot_SAA = sim.get_metabolite_time_series()["saa_e"].iloc[-1]
        tot_RA = sim.get_metabolite_time_series()["rosma_e"].iloc[-1]

        tot_BM_list.append(tot_BM)
        tot_CAs.append(tot_CA)
        tot_SAAs.append(tot_SAA)
        tot_RAs.append(tot_RA)
        inoc_ratio_list.append(inocculation_ratio)

    except Exception as e:
        print("not able to process:", inocculation_ratio)
        print(f"An exception occurred: {e}")

results_df = pd.DataFrame({'inoculation_ratio': inoc_ratio_list, 'total_biomass': tot_BM_list, 'total_RA': tot_RAs, 'total_CA': tot_CAs, 'total_SAA': tot_SAAs})

results_df.to_csv("../results/fig_3D", index=False)
