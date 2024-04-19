"""Script for running analysis to create "3D" production envelope for yield / BM / tot RA produced."""

import pandas as pd
import numpy as np
from tqdm import tqdm
from cobra.io import read_sbml_model
from modify_GEM import add_ratio_constraint_cobra
from dfba_comets import simulate_coculture

# ----- read cobrapy models -----

RAU2_cobra = read_sbml_model("../GEMs/RAU2.xml")
RAD4_cobra = read_sbml_model("../GEMs/RAD4.xml")

# ----- add yield ratio constraint -----

# knock out GLCtex_copy1 in both models, just so I have one less thing to worry about...
RAU2_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)
RAD4_cobra.reactions.get_by_id("GLCtex_copy1").bounds = (0.0, 0.0)

# ---- define yield values to iterate over ------
initial_yields = np.array([0.0304, 0.0362, 0.0128]) #from exp. values
yield_factors = np.arange(6.2, 12.2, 0.2)

# --- run analysis for each yield ---

biomass_list = []
RA_list = []
time_list = []
yield_factor_list = []

for yield_factor in tqdm(yield_factors):

    CA_yield, SAA_yield, RA_yield, = initial_yields * yield_factor

    # ----- add yield ratio constraint -----
    add_ratio_constraint_cobra(RAU2_cobra, "34DHCINMt", "GLCtex_copy2", CA_yield);
    add_ratio_constraint_cobra(RAD4_cobra, "SAAt", "GLCtex_copy2", SAA_yield);
    add_ratio_constraint_cobra(RAD4_cobra, "RAt", "GLCtex_copy2", RA_yield);

    try:
        sim = simulate_coculture(RAU2_cobra, RAD4_cobra, initial_pop=2.e-3, km_glc_adj=1000, glc_vmax_adjustment=0.45)

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
results_df.to_csv("../results/prod_env_yields_coculture.csv", index=False)