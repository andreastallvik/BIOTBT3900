"""
Result plots.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_relative_abundance_glc_xyl(results_df, cal11_ra, sal11_ra, mam3_ra):

    plot_df = results_df.explode(["CAL11", "SAL11", "MAM3"])
    df_melt = plot_df.melt('frac', var_name='strain', value_name='species abundance')

    sns.lineplot(data=df_melt, x="frac", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal11_ra, linestyle='--', color='blue')
    plt.axhline(y=sal11_ra, linestyle='--', color='orange')
    plt.axhline(y=mam3_ra, linestyle='--', color='green')

    plt.xlabel("fraction of optimal solution")


def plot_relative_abundance_RA_prod_glc_xyl(results_df, cal11_ra, sal11_ra, mam3_ra):
    
    plot_df = results_df.explode(["CAL11", "SAL11", "MAM3"])
    df_melt = plot_df.drop(columns=["RA_prod_rate"]).melt('RA_percentage', var_name='strain', value_name='species abundance')

    sns.lineplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal11_ra, linestyle='--', color='blue')
    plt.axhline(y=sal11_ra, linestyle='--', color='orange')
    plt.axhline(y=mam3_ra, linestyle='--', color='green')

    plt.xlabel("Percentage of maximal RA production")

    plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()


def plot_biomass_time_course(sim_results, exp_data):

    """ sim_results should be a dataframe outputted from sim.total_biomass """

    biomass_df = sim_results.copy()
    
    # convert from cycles to time, assuming a simulation rate of 0.1 hours
    biomass_df["time"] = biomass_df["cycle"]*0.1
    biomass_df.drop(columns=["cycle"], inplace=True)

    # add total biomass as a column
    biomass_df["total_BM"] = biomass_df["CAL11"] + biomass_df["SAL11"] + biomass_df["MAM3"]

    # prepare df for plotting
    plot_df = biomass_df.melt('time', var_name='variable', value_name='biomass')

    # plot
    sns.lineplot(data=plot_df, x="time", y="biomass", hue="variable")

    plt.scatter(x='time', y='biomass', data=exp_data, color = "red", marker=".", label="measured total biomass")
    plt.legend()
    plt.xlabel("time (h)")
    plt.ylabel("biomass (gDW)")


def plot_relative_abundance_time_course(sim_results, exp_data):

    """ sim_results should be a dataframe outputted from sim.total_biomass """
    
    biomass_df = sim_results.copy()

    total_BM = biomass_df["CAL11"] + biomass_df["SAL11"] + biomass_df["MAM3"]
    CAL11_frac = biomass_df["CAL11"] / total_BM
    SAL11_frac = biomass_df["SAL11"] / total_BM
    MAM3_frac = biomass_df["MAM3"] / total_BM

    df = pd.concat([total_BM, CAL11_frac, SAL11_frac, MAM3_frac], axis=1, keys=["total_BM", "CAL11", "SAL11", "MAM3"])
    df["time"] = df.index * 0.1
    plot_df = df.melt(id_vars="time", value_vars=["CAL11", "SAL11", "MAM3"], value_name="relative_abundance", var_name="strain")

    sns.lineplot(plot_df, x="time", y="relative_abundance", hue="strain", hue_order=["CAL11", "SAL11", "MAM3"])
    sns.scatterplot(data=exp_data, x='time', y='subpopulation_percentage', hue='strain', hue_order=["CAL11", "SAL11", "MAM3"], markers=".")
    
    plt.ylabel("Subpopulation fraction")
    plt.xlabel("Time (h)")


def plot_production_time_course(sim_results, exp_data):

    """ sim_results should be a dataframe outputted from sim.get_metabolite_time_series() """

    
    products_df = sim_results.copy()

    keep_cols = ["saa_e", "34dhcinm_e", "rosma_e", "cycle"]

    products_df = products_df[keep_cols].copy()
    products_df["time"] = products_df["cycle"]*0.1
    products_df.drop(columns=["cycle"], inplace=True)

    products_df.rename(columns={"saa_e": "SAA", "34dhcinm_e": "CA", "rosma_e": "RA"}, inplace=True)

    plot_df = products_df.melt(id_vars="time", value_vars=["SAA", "CA", "RA"], value_name="concentration", var_name="product")

    fig, ax1 = plt.subplots()

    # lineplot on the first y-axis (ax1)
    sns.lineplot(data=plot_df, x="time", y="concentration", hue="product", hue_order=["CA", "SAA", "RA"], ax=ax1)
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("simulated concentration (mmol/L)")

    # create a second y-axis (ax2) and plot scatterplot on it
    ax2 = ax1.twinx()
    sns.scatterplot(data=exp_data, x="time", y="mmol_per_L", hue="product", hue_order=["CA", "SAA", "RA"], ax=ax2, markers=".")
    ax2.set_ylabel("measured concentration (mmol/L)")

    plt.show()


def plot_inoculum_substrate():
    pass



def save_fig(filepath):
    """Save figure as vector image"""
    #TODO: implement
    pass