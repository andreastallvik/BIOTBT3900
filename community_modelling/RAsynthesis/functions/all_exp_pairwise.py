"""
Script for running 36 experiments with all inoculation ratio / glucose:xylose ratio combinations, 
in order to re-create fig. 6 from the RA paper.
Saves result as csv file in /results folder.
NOTE: Filepaths are hardocded, so scripts needs to be run from within the /function directory for anything to work.
"""

import itertools
#from dfba_comets import run_dfba
import pandas as pd
from tqdm import tqdm
import data_analysis
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_xyl_glc_triculture


## load and prepare models
print("loading exp. data")

# load exp. data yields, but magnified
relative_abundance_df = pd.read_csv("../exp_data/subpop_data_xyl_glc.csv")
products_df = pd.read_csv("../exp_data/conc_data_xyl_glc.csv")
od600_df = pd.read_csv("../exp_data/od600_xyl_glc.csv")
subpop_df, conc_df, OD_df = data_analysis.process_data(relative_abundance_df, products_df, od600_df)
CA_yield, SAA_yield, RA_yield = data_analysis.get_yields_glc_xyl(conc_df)
# CA_yield = CA_yield * 10
# SAA_yield = SAA_yield * 10
# RA_yield = RA_yield * 10

print("loading and preparing GEMs")

# ----- read cobrapy models -----

CAL11_cobra = read_sbml_model("../GEMs/CAL11.xml")
SAL11_cobra = read_sbml_model("../GEMs/SAL11.xml")
MAM3_cobra = read_sbml_model("../GEMs/MAM3.xml")

# ----- add yield ratio constraint -----

# knock out GLCtex_copy1 in SAL11, just so I have one less thing to worry about...
SAL11_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

add_ratio_constraint_cobra(CAL11_cobra, "34DHCINMt", "XYLtex", CA_yield);
add_ratio_constraint_cobra(SAL11_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
add_ratio_constraint_cobra(MAM3_cobra, "RAt", "XYLtex", RA_yield);


inocculation_ratios = [
    (2,3,1),
    (1,3,1),
    (3,3,1),
    (1,1,1),
    (1,2,1),
    (2,2,1),
    (2,1,1),
    (3,1,1),
    (3,2,1)
]

glucose_xylose_ratios = [
    (1, 4),
    (2, 3),
    (3, 2),
    (4, 1)
]

def convert_to_mmol(glc_xyl_ratio):

    MM_glc = 180.156
    MM_xyl = 150.13

    glc_g_L = glc_xyl_ratio[0]
    xyl_g_L = glc_xyl_ratio[1]

    glc_g = glc_g_L * 0.1
    xyl_g = xyl_g_L * 0.1

    glc_m = glc_g / MM_glc
    xyl_m = xyl_g / MM_xyl

    glc_mmol = glc_m * 1000
    xyl_mmol = xyl_m * 1000

    return round(glc_mmol, 2), round(xyl_mmol, 2)

print("setting up experiment")

all_experiment_combinations = list(itertools.product(inocculation_ratios, glucose_xylose_ratios))


tot_BM_list = []
tot_RAs = []
inoc_ratio_list = []
glc_xyl_ratio_list = []

for experiment in tqdm(all_experiment_combinations):
    inocculation_ratio = experiment[0]
    glucose_xylose_ratio = convert_to_mmol(experiment[1])

    try:
        #sim = run_dfba(initial_pop_ratio=inocculation_ratio, glc_xyl_mmol=glucose_xylose_ratio, RA_lb=0.5)

        sim = simulate_xyl_glc_triculture(CAL11_cobra, SAL11_cobra, MAM3_cobra, 
                                          initial_pop=3.9e-3, initial_pop_ratio=inocculation_ratio, 
                                          glc_xyl_mmol=glucose_xylose_ratio, adjust_atp_requirements=True)

        tot_BM = sum(sim.total_biomass.drop(columns=["cycle"], inplace=False).iloc[-1])
        tot_RA = sim.get_metabolite_time_series()["rosma_e"].iloc[-1]

        tot_BM_list.append(tot_BM)
        tot_RAs.append(tot_RA)
        inoc_ratio_list.append(inocculation_ratio)
        glc_xyl_ratio_list.append(experiment[1])

    except Exception as e:
        print("not able to process combination:", experiment)
        print(f"An exception occurred: {e}")

results_df = pd.DataFrame({'inoculation_ratio': inoc_ratio_list, 'glc_xyl_ratio': glc_xyl_ratio_list, 'total_biomass': tot_BM_list, 'total_RA': tot_RAs})

results_df.to_csv("../results/fig_6_new_params", index=False)
