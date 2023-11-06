"""Script for running analysis to create "3D" production envelope for yield / BM / tot RA produced."""

import pandas as pd
import numpy as np
from tqdm import tqdm
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_glc_triculture

# ----- read cobrapy models -----

CAL2_cobra = read_sbml_model("../GEMs/CAL2.xml")
SAL9_cobra = read_sbml_model("../GEMs/SAL9.xml")
MAM2_cobra = read_sbml_model("../GEMs/MAM2.xml")

# knock out GLCtex_copy1 in all organisms, just so I have one less thing to worry about when flux coupling
CAL2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
SAL9_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
MAM2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

# ---- define yield values to iterate over ------
initial_yields = np.array([0.0448, 0.115, 0.0600]) #from exp. values
yield_factors = np.arange(2.5, 7.6, 0.2)

# --- run analysis for each yield ---

biomass_list = []
RA_list = []
time_list = []
yield_factor_list = []

for yield_factor in tqdm(yield_factors):

    CA_yield, SAA_yield, RA_yield, = initial_yields * yield_factor

    # ----- add yield ratio constraint -----
    add_ratio_constraint_cobra(CAL2_cobra, "34DHCINMt", "GLCtex_copy2", CA_yield)
    add_ratio_constraint_cobra(SAL9_cobra, "SAAt", "GLCtex_copy2", SAA_yield)
    add_ratio_constraint_cobra(MAM2_cobra, "RAt", "GLCtex_copy2", RA_yield)

    try:
        sim = simulate_glc_triculture(CAL2_cobra, SAL9_cobra, MAM2_cobra, initial_pop=2.7e-3, initial_pop_ratio=(2, 3, 1), 
                              adjust_atp_requirements=True, km_glc_adj=1000, glc_vmax_adjustment=0.45)

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
        break


biomass_df = pd.concat(biomass_list, axis=0)
RA_df = pd.concat(RA_list, axis=0)
time_df = pd.concat(time_list, axis=0)
yield_factor_df = pd.concat(yield_factor_list, axis=0)

results_df = pd.concat([biomass_df, RA_df, time_df, yield_factor_df], axis=1, ignore_index=True, keys=["biomass", "RA", "time", "yield_factor"])
results_df = results_df.rename({"0": "biomass", "1": "RA", "2": "time","3": "yield_factor"}, axis=1)
results_df.to_csv("../results/prod_env_yields_glc.csv", index=False)