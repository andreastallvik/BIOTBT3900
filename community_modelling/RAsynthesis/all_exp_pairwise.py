
import itertools
from dfba_comets import run_dfba
import pandas as pd
from tqdm import tqdm

inocculation_ratios = [
    (2,3,1),
    (1,3,1),
    (3,3,1),
    (1,1,1),
    (1,2,1),
    (2,2,1),
    (2,1,1), #
    (3,1,1), #
    (3,2,1) #
]

glucose_xylose_ratios = [
    (1, 4),
    (2, 3),
    (3, 2), #
    (4, 1) #
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


all_experiment_combinations = list(itertools.product(inocculation_ratios, glucose_xylose_ratios))


tot_BM_list = []
tot_RAs = []
inoc_ratio_list = []
glc_xyl_ratio_list = []

for experiment in tqdm(all_experiment_combinations):
    inocculation_ratio = experiment[0]
    glucose_xylose_ratio = convert_to_mmol(experiment[1])

    try:
        sim = run_dfba(initial_pop_ratio=inocculation_ratio, glc_xyl_mmol=glucose_xylose_ratio)

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

results_df.to_csv("all_exp_results_full", index=False)
