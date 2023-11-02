"""
Functions for processing the experimental data, adding useful columns, and printing uselful stats.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def process_data(subpop_df, conc_df, OD_df):

    # remove negative values
    subpop_df = remove_neg_values(subpop_df, "subpopulation_percentage")
    conc_df = remove_neg_values(conc_df, "concentration")
    OD_df = remove_neg_values(OD_df, "time")

    # add row for converted OD meaurements
    OD_df["biomass"] = OD600_to_BM(OD_df["OD600"])

    # add row for converted mmol for subpop data
    conc_df['mmol'] = conc_df.apply(mg_per_L_to_mmol, axis=1)
    conc_df['mmol_per_L'] = conc_df.apply(mg_per_L_to_mmol_per_L, axis=1)

    return subpop_df, conc_df, OD_df


def remove_neg_values(df, var):
    """negative values are small errors from the data extrapolation from the figure - reset these to 0"""
    result_df = df.copy()
    result_df[var] = result_df[var].apply(lambda x: max(0, x))
    return result_df


def OD600_to_BM(OD600_measurement):
    """Conversion of OD600 to biomass for E. coli. source: https://doi.org/10.1016/j.ymben.2016.05.006"""
    BM_per_L = 0.31*OD600_measurement
    BM = BM_per_L*0.1
    return BM


def mg_per_L_to_mmol(row):

    # molar masses taken from KEGG
    MM_CA = 180.1574
    MM_SAA = 198.1727
    MM_RA = 360.3148

    mg_per_L = row['concentration']
    product = row['product']

    if product == 'CA':
        molar_mass = MM_CA
    elif product == 'SAA':
        molar_mass = MM_SAA
    elif product == 'RA':
        molar_mass = MM_RA

    mg = mg_per_L * 0.1
    mmol = mg / molar_mass
    return mmol


def mg_per_L_to_mmol_per_L(row):

    # molar masses taken from KEGG
    MM_CA = 180.1574
    MM_SAA = 198.1727
    MM_RA = 360.3148

    mg_per_L = row['concentration']
    product = row['product']

    if product == 'CA':
        molar_mass = MM_CA
    elif product == 'SAA':
        molar_mass = MM_SAA
    elif product == 'RA':
        molar_mass = MM_RA

    mmol_per_L = mg_per_L / molar_mass
    return mmol_per_L


def get_yields_glc_xyl(products_df):

    # get final concentrations (in mg/L)
    CA_c = products_df[products_df["product"] == "CA"]["concentration"].iloc[-1]
    SAA_c = products_df[products_df["product"] == "SAA"]["concentration"].iloc[-1]
    RA_c = products_df[products_df["product"] == "RA"]["concentration"].iloc[-1]

    # molar masses taken from KEGG
    MM_CA = 180.1574
    MM_SAA = 198.1727
    MM_RA = 360.3148

    # convert to mmol
    CA_mmol = CA_c * 0.1 / MM_CA
    SAA_mmol = SAA_c * 0.1 / MM_SAA
    RA_mmol = RA_c * 0.1 / MM_RA

    # calculate the TOTAL CA and SAA amount produced (since they are converted to RA, there is more produced than is accumelated at exp. end)
    tot_CA_mmol = CA_mmol + RA_mmol
    tot_SAA_mmol = SAA_mmol + RA_mmol

    # covert to g/L to get the same units as substrate | sequentially: mmol -> mol -> g -> g/L
    tot_CA_c = tot_CA_mmol * 0.001 * MM_CA / 0.1
    tot_SAA_c = tot_SAA_mmol * 0.001 * MM_SAA / 0.1
    tot_RA_c = RA_mmol * 0.001 * MM_RA / 0.1

    # # convert mg/L to g/L (same units as substrate)
    # CA_c = CA_c * 0.0001
    # SAA_c = SAA_c * 0.0001
    # RA_c = RA_c * 0.0001

    #NOTE: assuming that CA and RA module consume glucose in proportion with their steady-state relative abundance
    CA_yield = tot_CA_c / 1.23
    SAA_yield = tot_SAA_c / 3
    RA_yield = tot_RA_c / 0.77

    print("CA yield", CA_yield, "g CA per g xylose")
    print("SAA yield", SAA_yield, "g SAA per g glucose")
    print("RA yield", RA_yield, "g RA per g xylose")

    return CA_yield, SAA_yield, RA_yield


def get_yields_glc(products_df):

    # get final concentrations (in mg/L)
    CA_c = products_df[products_df["product"] == "CA"]["concentration"].iloc[-1]
    SAA_c = products_df[products_df["product"] == "SAA"]["concentration"].iloc[-1]
    RA_c = products_df[products_df["product"] == "RA"]["concentration"].iloc[-1]

    # molar masses taken from KEGG
    MM_CA = 180.1574
    MM_SAA = 198.1727
    MM_RA = 360.3148

    # convert to mmol
    CA_mmol = CA_c * 0.1 / MM_CA
    SAA_mmol = SAA_c * 0.1 / MM_SAA
    RA_mmol = RA_c * 0.1 / MM_RA

    # calculate the TOTAL CA and SAA amount produced (since they are converted to RA, there is more produced than is accumelated at exp. end)
    tot_CA_mmol = CA_mmol + RA_mmol
    tot_SAA_mmol = SAA_mmol + RA_mmol

    # covert to g/L to get the same units as substrate | sequentially: mmol -> mol -> g -> g/L
    tot_CA_c = tot_CA_mmol * 0.001 * MM_CA / 0.1
    tot_SAA_c = tot_SAA_mmol * 0.001 * MM_SAA / 0.1
    tot_RA_c = RA_mmol * 0.001 * MM_RA / 0.1

    # NOTE: assuming glucose in consumed in proportion with steady-state relative abundance
    CA_yield = tot_CA_c / 2.5
    SAA_yield = tot_SAA_c / 0.8
    RA_yield = tot_RA_c / 1.7

    print("CA yield", CA_yield, "g CA per g glucose")
    print("SAA yield", SAA_yield, "g SAA per g glucose")
    print("RA yield", RA_yield, "g RA per g glucose")

    return CA_yield, SAA_yield, RA_yield


def get_yields_coculture(CA_c = 60, SAA_c = 73, RA_c = 32):

    # molar masses taken from KEGG
    MM_CA = 180.1574
    MM_SAA = 198.1727
    MM_RA = 360.3148

    # convert to mmol
    CA_mmol = CA_c * 0.1 / MM_CA
    SAA_mmol = SAA_c * 0.1 / MM_SAA
    RA_mmol = RA_c * 0.1 / MM_RA

    # calculate the TOTAL CA and SAA amount produced (since they are converted to RA, there is more produced than is accumelated at exp. end)
    tot_CA_mmol = CA_mmol + RA_mmol
    tot_SAA_mmol = SAA_mmol + RA_mmol

    # covert to g/L to get the same units as substrate | sequentially: mmol -> mol -> g -> g/L
    tot_CA_c = tot_CA_mmol * 0.001 * MM_CA / 0.1
    tot_SAA_c = tot_SAA_mmol * 0.001 * MM_SAA / 0.1
    tot_RA_c = RA_mmol * 0.001 * MM_RA / 0.1

    #NOTE: making the (likely incorrect) assumption that each strain consumes 2.5 g/L glucose each
    CA_yield = tot_CA_c / 2.5
    SAA_yield = tot_SAA_c / 2.5
    RA_yield = tot_RA_c / 2.5

    print("CA yield", CA_yield, "g CA per g glucose")
    print("SAA yield", SAA_yield, "g SAA per g glucose")
    print("RA yield", RA_yield, "g RA per g glucose")

    return CA_yield, SAA_yield, RA_yield


def print_production_stats_triculture(products_df):
    print("Final RA concentration (mg/L):", products_df[products_df["product"] == "RA"]["concentration"].iloc[-1])
    print("Final RA amount (mmol):", products_df[products_df["product"] == "RA"]["mmol"].iloc[-1])

    print("Final SAA concentration (mg/L):", products_df[products_df["product"] == "SAA"]["concentration"].iloc[-1])
    print("Final SAA amount (mmol):", products_df[products_df["product"] == "SAA"]["mmol"].iloc[-1])

    print("Final CA concentration (mg/L):", products_df[products_df["product"] == "CA"]["concentration"].iloc[-1])
    print("Final CA amount (mmol):", products_df[products_df["product"] == "CA"]["mmol"].iloc[-1])


def get_relative_abundance_stats_triculture(subpop_df):
    """Get relative abundance at steady-state, and print out. Returns relative abundance in order CAL, SAL, MAM."""

    sal11_ra = subpop_df[subpop_df["strain"] == "SAL11"]["subpopulation_percentage"].iloc[-1]
    cal11_ra = subpop_df[subpop_df["strain"] == "CAL11"]["subpopulation_percentage"].iloc[-1]
    mam3_ra = subpop_df[subpop_df["strain"] == "MAM3"]["subpopulation_percentage"].iloc[-1]

    print("Relative SAl11 abundance at steady-state:", sal11_ra)
    print("Relative CAL11 abundance at steady-state:", cal11_ra)
    print("Relative MAM3 abundance at steady-state:", mam3_ra)    

    return cal11_ra, sal11_ra, mam3_ra


def get_relative_abundance_stats_triculture_glc(subpop_df):
    """Get relative abundance at steady-state, and print out. Returns relative abundance in order CAL, SAL, MAM. For CAL2:SAL9:MAM2."""

    sal9_ra = subpop_df[subpop_df["strain"] == "SAL9"]["subpopulation_percentage"].iloc[-1]
    cal2_ra = subpop_df[subpop_df["strain"] == "CAL2"]["subpopulation_percentage"].iloc[-1]
    mam2_ra = subpop_df[subpop_df["strain"] == "MAM2"]["subpopulation_percentage"].iloc[-1]

    print("Relative SAl9 abundance at steady-state:", sal9_ra)
    print("Relative CAL2 abundance at steady-state:", cal2_ra)
    print("Relative MAM2 abundance at steady-state:", mam2_ra)    

    return cal2_ra, sal9_ra, mam2_ra


def get_growth_curves(od_df, subpop_df, strains):
    """Calculates individual strain biomass measurements using total OD and relative abundance. 
    Assumes "time" and "strain" columns in inputs are ordered (accending time, clustered by strain)."""
    
    strain_bm = []

    for strain in strains:
        bm = od_df["biomass"].values * subpop_df[subpop_df["strain"] == strain]["subpopulation_percentage"].values
        strain_bm.append(bm.transpose())

    bm_df = pd.DataFrame(np.array(strain_bm).T, columns=strains)
    bm_df["time"] = od_df["time"]

    return bm_df.melt(id_vars="time", value_name="biomass", var_name="strain")