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
#yield_factors = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5])
#yield_factors = np.arange(1.3, 3.3, 0.2)
#yield_factors = np.arange(1.4, 8.0, 0.4) # stopped being possible at 5.4
yield_factors = np.arange(2.5, 5.5, 0.1)

# --- run analysis for each yield ---

biomass_list = []
RA_list = []
time_list = []
yield_factor_list = []

for yield_factor in tqdm(yield_factors):

    CA_yield, SAA_yield, RA_yield, = initial_yields * yield_factor

    # ----- add yield ratio constraint -----
    add_ratio_constraint_cobra(CAL11_cobra, "34DHCINMt", "XYLtex", CA_yield);
    add_ratio_constraint_cobra(SAL11_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
    add_ratio_constraint_cobra(MAM3_cobra, "RAt", "XYLtex", RA_yield);

    try:
        sim = simulate_xyl_glc_triculture(CAL11_cobra, SAL11_cobra, MAM3_cobra, 
                                          initial_pop=3.9e-3, adjust_atp_requirements=True)

        biomass = sim.total_biomass.drop(columns=["cycle"], inplace=False).sum(axis="columns")
        time = sim.total_biomass["cycle"] * 0.1
        RA = sim.get_metabolite_time_series()["rosma_e"]

        biomass_list.append(biomass)
        time_list.append(time)
        RA_list.append(RA)
        yield_factor_list.append(pd.Series([yield_factor] * len(time)))

    except Exception as e:
        print("not able to process yield factor:", yield_factor)
        print(f"An exception occurred: {e}")


biomass_df = pd.concat(biomass_list, axis=0)
RA_df = pd.concat(RA_list, axis=0)
time_df = pd.concat(time_list, axis=0)
yield_factor_df = pd.concat(yield_factor_list, axis=0)

results_df = pd.concat([biomass_df, RA_df, time_df, yield_factor_df], axis=1, ignore_index=True, keys=["biomass", "RA", "time", "yield_factor"])
results_df = results_df.rename({"0": "biomass", "1": "RA", "2": "time","3": "yield_factor"}, axis=1)
results_df.to_csv("../results/prod_env_yields_tri_xyl_glc_new_params.csv", index=False)