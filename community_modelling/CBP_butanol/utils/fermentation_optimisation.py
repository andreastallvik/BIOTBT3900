"""
Script to search the parameter space for the optimal fermentation conditions of:
inoculation time, inoculation ratio, and xylan concentration, to get highest butanol production.

NOTE: file-paths are absolute, so the script needs to be run from the CBP_butanol directory.
"""

from comets_functions import sequental_com, collapse_sequential_sim
import pandas as pd
from cobra.io import read_sbml_model
from flux_coupling import add_ratio_constraint_cobra
from kinetic_params import KINETIC_PARAMS
import datetime
from tqdm import tqdm

# constants
VOLUME = 0.05
MM_XYLAN8 = 1201.04
timestamp = datetime.datetime.now().strftime("%b_%d_%H%M")
FILEPATH = f"grid_search_results/grid_search_result_{timestamp}.csv"

# ---------------- load models ----------------
print("loading models...")
nj4 = read_sbml_model("GEMs/NJ4_curated.xml")
m5 = read_sbml_model("GEMs/M5_curated.xml")

# ---------------- add constraints to models ----------------
print("adding constraints...")

# make the reactions in the ABE pathway irreversible

reactions = ["POR_syn",
            "ACACT1r",
            "HACD1",
            "ECOAH1",
            "ACOAD1fr",
            "ACOAD1",
            "BTCOARx",
            "PBUTT",
            "ADCi",
            "PTAr"]

reverse_reactions = ["ALCD4", "BUTKr", "BUTCT2", "ACKr", "ACACCT", "ACALD"]

for rx in reactions:
    nj4.reactions.get_by_id(rx).bounds = (0, 1000)

for rx in reverse_reactions:
    nj4.reactions.get_by_id(rx).bounds = (-1000, 0)

# knock out reactions for xylan uptake

xylan_rx = ["XYLANabc", "GLCURS1"]

for rx in xylan_rx:
    nj4.reactions.get_by_id(rx).bounds = (0, 0)

# knock out reactions for acetate and butyrate production
nj4.reactions.ACtr.bounds = (0, 1000)
nj4.reactions.BUTt.bounds = (0, 1000)

# ac / but flux uptake flux couplued to adhere to exp. measurement
add_ratio_constraint_cobra(nj4, "BUTt" , "ACtr",  0.42, r_num_reverse=False, r_den_reverse=False)

# acetone / butanol production flux coupling, constrained to exp. measurement
add_ratio_constraint_cobra(nj4, "BTOHt" , "ACEt",  2.68, r_num_reverse=True, r_den_reverse=False)

# knock out reactions for xylose uptake
uptake_KO = ["XYLANabc", "XYLabc", "XYLtex"]
for rx in uptake_KO:
    m5.reactions.get_by_id(rx).bounds = (0, 0)

# restrict the max rate of xylose uptake
m5.reactions.XYLt2.bounds = (0, 0.25)

# constrain butyrate / acetate production ratio to exp. measurement
add_ratio_constraint_cobra(m5, "BUTt" , "ACtr",  0.71, r_num_reverse=False, r_den_reverse=False)

# ---------------- get the base medium dict ----------------
media_db = pd.read_csv("medium.tsv", sep="\t")
m5_med = media_db[media_db["medium"] == "m5_med"]

UNLIMITED_METABOLITES = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'k_e', 'h2o_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'pi_e', 'so4_e', 'zn2_e']

metabolite_list = [str(m+"_e") for m in m5_med["compound"].tolist()]
limited_metabolites = set(metabolite_list) - set(UNLIMITED_METABOLITES)
medium = {k:0.5 for k in limited_metabolites}

# ---------------- grid search ----------------
print("running the grid search...")
inoculation_times = [60, 70]
inoculation_ratios = [1, 1.5]
xylan_concentrations = [60, 70]

results = pd.DataFrame(columns=['inoc_time', 'inoc_ratio', 'xylan_conc', 'butanol'])

# set up progress bar
total = len(inoculation_times) * len(inoculation_ratios) * len(xylan_concentrations)
pbar = tqdm(total=total)

for inoculation_time in inoculation_times:
    for inoculation_ratio in inoculation_ratios:
        for xylan_concentration in xylan_concentrations:

            # calculate the xylan amount in mmol
            xylan = (xylan_concentration * VOLUME / MM_XYLAN8) * 1000
            # update the medium
            medium["xylan_e"] = xylan

            try:
                # run the simulation
                first_sim, second_sim = sequental_com(m5, nj4, init_medium=medium, kinetic_params=KINETIC_PARAMS, inoc_time=inoculation_time, inoc_ratio=inoculation_ratio)
                
                # collapse the results
                bm, met, fluxes = collapse_sequential_sim(first_sim, second_sim)

                # get final butanol titer
                if "btoh_e" in met.columns:
                    btoh = met["btoh_e"].iloc[-1]
                else:
                    btoh = 0

            except Exception as e:
                print("Exception for inoculation_time: " + str(inoculation_time) + ", inoculation_ratio: " + str(inoculation_ratio) + ", xylan_concentration: " + str(xylan_concentration))
                print(e)
                btoh = float("NaN")

            # save the results
            current_res = pd.DataFrame({'inoc_time': inoculation_time, 'inoc_ratio': inoculation_ratio, 'xylan_conc': xylan_concentration, 'butanol': btoh} , index=[0])
            results = pd.concat([results, current_res], ignore_index=True)

            # write the dataframe to file as csv, doing this so that some results are available even if program is quit prematurely
            results.to_csv(FILEPATH, index=False)

            # update the progress bar
            pbar.update()

# close the progress bar
pbar.close()