"""Script for running analysis to create "3D" production envelope for yield / BM / tot RA produced."""

import pandas as pd
import numpy as np
from tqdm import tqdm
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_xyl_glc_triculture

# ----- read cobrapy models -----

CAL11_cobra = read_sbml_model("../GEMs/CAL11.xml")
SAL11_cobra = read_sbml_model("../GEMs/SAL11.xml")
MAM3_cobra = read_sbml_model("../GEMs/MAM3.xml")

# knock out GLCtex_copy1 in SAL11
SAL11_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

# ---- define yield values to iterate over ------
initial_yields = np.array([0.118, 0.0431, 0.223]) #from exp. values
yield_factors = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5])


# --- run analysis for each yield ---

tot_BM_list = []
tot_RAs = []
tot_yield_factors = []

for yield_factor in tqdm(yield_factors):

    CA_yield, SAA_yield, RA_yield, = initial_yields * yield_factor

    # ----- add yield ratio constraint -----
    add_ratio_constraint_cobra(CAL11_cobra, "34DHCINMt", "XYLtex", CA_yield);
    add_ratio_constraint_cobra(SAL11_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
    add_ratio_constraint_cobra(MAM3_cobra, "RAt", "XYLtex", RA_yield);

    try:
        sim = simulate_xyl_glc_triculture(CAL11_cobra, SAL11_cobra, MAM3_cobra, 
                                          initial_pop=2.e-3, adjust_atp_requirements=True)

        tot_BM = sum(sim.total_biomass.drop(columns=["cycle"], inplace=False).iloc[-1])
        tot_RA = sim.get_metabolite_time_series()["rosma_e"].iloc[-1]

        tot_BM_list.append(tot_BM)
        tot_RAs.append(tot_RA)
        tot_yield_factors.append(yield_factor)

    except Exception as e:
        print("not able to process yield factor:", yield_factor)
        print(f"An exception occurred: {e}")

results_df = pd.DataFrame({'yield_factor': tot_yield_factors, 'total_biomass': tot_BM_list, 'total_RA': tot_RAs})

results_df.to_csv("../results/prod_env_yields_tri_xyl_glc.csv", index=False)