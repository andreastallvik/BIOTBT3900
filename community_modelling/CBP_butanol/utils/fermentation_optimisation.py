"""
Script to search the parameter space for the optimal fermentation conditions of:
inoculation time, inoculation ratio, and xylan concentration, to get highest butanol production.

NOTE: file-paths are absolute, so the script needs to be run from the CBP_butanol directory.
"""

from comets_functions import sequential_with_switch, collapse_three_sim
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

# parameter grid
inoculation_times = [24, 48, 72]
inoculation_ratios = [0.5, 1, 1.33]
xylan_concentrations = [40, 60, 80]

# ---------------- load models ----------------
print("loading models...")
nj4 = read_sbml_model("GEMs/NJ4_curated.xml")
m5 = read_sbml_model("GEMs/M5_curated.xml")

# ---------------- add constraints to models ----------------
print("adding constraints...")

# M5:

uptake_KO = ["XYLANabc", "XYLabc", "XYLtex"]

for rx in uptake_KO:
    m5.reactions.get_by_id(rx).bounds = (0, 0)

# restrict rate of xylose uptake
m5.reactions.XYLt2.bounds = (0, 0.4)

# restrict uptake of butanol
m5.reactions.BTOHt.bounds = (-1000, 0)

# flux coupling constraint forcing but/ac production to exp. values
add_ratio_constraint_cobra(m5, "BUTt" , "ACtr",  0.71, r_num_reverse=False, r_den_reverse=False)

# nj4
# define the Specific Proton Flux (SPF) property

h_membrane_rx = [r.id for r in nj4.metabolites.h_e.reactions if "EX" not in r.id]

neg_stoich = []
pos_stoich = []

for rx in h_membrane_rx:
    stociometry = {met.id:coeff for met, coeff in nj4.reactions.get_by_id(rx).metabolites.items()}    
    if stociometry["h_e"] < 0:
        pos_stoich.append(rx)
    elif stociometry["h_e"] > 0:
        neg_stoich.append(rx)

# as an objective
SPF_obj = nj4.problem.Objective(
        sum([nj4.reactions.get_by_id(rx).flux_expression for rx in pos_stoich]) - sum([nj4.reactions.get_by_id(rx).flux_expression for rx in neg_stoich]),
        direction="max")

nj4_acido = nj4.copy()
nj4_solvento = nj4.copy()

# restrict reaction reversibility

reactions = ["ACACT1r", "HACD1", "ECOAH1", "ACOAD1fr", "ACOAD1", "BTCOARx", "PBUTT", 
             "ADCi", "PTAr", "POR_syn", "FNOR", "FNRR","T2ECR", "BNOCA", #ABE pathway
             ]

reverse_reactions = ["ALCD4", "BUTKr", "BUTCT2", "ACKr", "ACACCT", "ACALD",
                      "HYDA", "HACD1i", "ACOAD1fr", "ACOAD1f", #ABE pathway
                      ]

KO_rx = ["XYLANabc", "GLCURS1"] #xylan uptake reactions

for rx in reactions:
    nj4_acido.reactions.get_by_id(rx).bounds = (0, 1000)
    nj4_solvento.reactions.get_by_id(rx).bounds = (0, 1000)

for rx in reverse_reactions:
    nj4_acido.reactions.get_by_id(rx).bounds = (-1000, 0)
    nj4_solvento.reactions.get_by_id(rx).bounds = (-1000, 0)

for rx in KO_rx:
    nj4_acido.reactions.get_by_id(rx).bounds = (0, 0)
    nj4_solvento.reactions.get_by_id(rx).bounds = (0, 0)

# TCA cycle
nj4_acido.reactions.SUCD2.bounds = (-1000, 1000)

# "closing off" reductive TCA for solventogenesis
nj4_solvento.reactions.SUCD2.bounds = (-1000, 0)
nj4_solvento.reactions.MDH.bounds = (-1000, 0)
# nj4_solvento.reactions.FUM.bounds = (-1000, 0)

# knock out reactions for acetate and butyrate production
nj4_solvento.reactions.ACtr.bounds = (0, 1000)
nj4_solvento.reactions.BUTt.bounds = (0, 1000)

# flux coupling but/ac production for acidogeneis to exp. value
add_ratio_constraint_cobra(nj4_acido, "BUTt" , "ACtr",  1.02, r_num_reverse=True, r_den_reverse=True)

# flux coupling btoh/acetone production for solventogenesis to exp. value
add_ratio_constraint_cobra(nj4_solvento, "BTOHt" , "ACEt",  2.68, r_num_reverse=True, r_den_reverse=False)

# add SPF objectove to solventogen model
nj4_solvento.objective = SPF_obj

# update some kinetic params
KINETIC_PARAMS["M5"]["km"]["EX_xylan8_e"] = 0.7
KINETIC_PARAMS["M5"]["km"]["EX_xyl__D_e"] = 10
KINETIC_PARAMS["NJ4"]["km"]["EX_xyl__D_e"] = 1

AA_uptake_rx = ["EX_val__L_e", "EX_arg__L_e", "EX_asp__L_e", "EX_dhptd_e", "EX_glu__L_e", "EX_ile__L_e", "EX_ser__L_e", 
                "EX_thr__L_e", "EX_ala__L_e", "EX_cys__L_e", "EX_gly_e", "EX_his__L_e", "EX_leu__L_e", 
                "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_tyr__L_e", "EX_trp__L_e", "EX_lys__L_e"]

for strain in KINETIC_PARAMS.keys():
    for rx in AA_uptake_rx:
        KINETIC_PARAMS[strain]["km"][rx] = 1

# ---------------- get the base medium dict ----------------
media_db = pd.read_csv("medium.tsv", sep="\t")
m5_med = media_db[media_db["medium"] == "m5_med"]

UNLIMITED_METABOLITES = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'k_e', 'h2o_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'pi_e', 'so4_e', 'zn2_e']

metabolite_list = [str(m+"_e") for m in m5_med["compound"].tolist()]
limited_metabolites = set(metabolite_list) - set(UNLIMITED_METABOLITES)
medium = {k:0.5 for k in limited_metabolites}
medium["xylan4_e"] = 0

# ---------------- grid search ----------------
print("running the grid search...")

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
                switch_time = inoculation_time + 24 #trying this
                first_sim, second_sim, third_sim = sequential_with_switch(m5=m5, nj4_acido=nj4_acido, nj4_solvento=nj4_solvento, init_medium=medium, 
                                                  kinetic_params=KINETIC_PARAMS, inoc_time=inoculation_time, inoc_ratio=inoculation_ratio, switch_time=switch_time)
                
                # collapse the results
                bm, met, fluxes = collapse_three_sim(first_sim, second_sim, third_sim)

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