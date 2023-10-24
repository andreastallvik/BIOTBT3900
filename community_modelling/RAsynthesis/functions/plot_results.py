"""
Result plots.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_relative_abundance_glc_xyl(results_df, cal11_ra, sal11_ra, mam3_ra, scatter=False):
    """Make relative abundance plot for CAL11:SAL11:MAM11 tri-culture"""

    plot_df = results_df.explode(["CAL11", "SAL11", "MAM3"])
    df_melt = plot_df.melt('frac', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="frac", y="species abundance", hue="strain")
    else:
        sns.lineplot(data=df_melt, x="frac", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal11_ra, linestyle='--', color='blue')
    plt.axhline(y=sal11_ra, linestyle='--', color='orange')
    plt.axhline(y=mam3_ra, linestyle='--', color='green')

    plt.xlabel("fraction of optimal solution")


def plot_relative_abundance_glc(results_df, cal2_ra, sal9_ra, mam2_ra, scatter=False):
    """Make relative abundance plot for CAL2:SAL9:MAM2 tri-culture"""

    plot_df = results_df.explode(["CAL2", "SAL9", "MAM2"])
    df_melt = plot_df.melt('frac', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="frac", y="species abundance", hue="strain")
    else:
        sns.lineplot(data=df_melt, x="frac", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal2_ra, linestyle='--', color='blue')
    plt.axhline(y=sal9_ra, linestyle='--', color='orange')
    plt.axhline(y=mam2_ra, linestyle='--', color='green')

    plt.xlabel("fraction of optimal solution")


def plot_relative_abundance_coculture(results_df, scatter=False):
    """Make relative abundance plot for CAL2:SAL9:MAM2 tri-culture"""

    plot_df = results_df.explode(["RAU2", "RAD4"])
    df_melt = plot_df.melt('frac', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="frac", y="species abundance", hue="strain")
    else:
        sns.lineplot(data=df_melt, x="frac", y="species abundance", hue="strain",orient="y")

    plt.xlabel("fraction of optimal solution")


def plot_relative_abundance_RA_prod_glc_xyl(results_df, cal11_ra, sal11_ra, mam3_ra, scatter=False):
    """Make relative abundance plot different production-rates of RA, for CAL11:SAL11:MAM3 tri-culture"""
    
    plot_df = results_df.explode(["CAL11", "SAL11", "MAM3"])
    df_melt = plot_df.drop(columns=["RA_prod_rate"]).melt('RA_percentage', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain")    
    else:
        sns.lineplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal11_ra, linestyle='--', color='blue')
    plt.axhline(y=sal11_ra, linestyle='--', color='orange')
    plt.axhline(y=mam3_ra, linestyle='--', color='green')

    plt.xlabel("Percentage of maximal RA production")

    plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()


def plot_relative_abundance_RA_prod_glc(results_df, cal2_ra, sal9_ra, mam2_ra, scatter=False):
    """Make relative abundance plot different production-rates of RA, for CAL2:SAL9:MAM2 tri-culture"""
    
    plot_df = results_df.explode(["CAL2", "SAL9", "MAM2"])
    df_melt = plot_df.drop(columns=["RA_prod_rate"]).melt('RA_percentage', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain")    
    else:
        sns.lineplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain",orient="y")

    #stipled lines with exp. steady-state species abundance
    plt.axhline(y=cal2_ra, linestyle='--', color='blue')
    plt.axhline(y=sal9_ra, linestyle='--', color='orange')
    plt.axhline(y=mam2_ra, linestyle='--', color='green')

    plt.xlabel("Percentage of maximal RA production")

    plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()


def plot_relative_abundance_RA_prod_coculture(results_df, scatter=False):
    """Make relative abundance plot different production-rates of RA, for RAU2:RAD4 co-culture"""
    
    plot_df = results_df.explode(["RAU2", "RAD4"])
    df_melt = plot_df.drop(columns=["RA_prod_rate"]).melt('RA_percentage', var_name='strain', value_name='species abundance')

    if scatter:
        sns.scatterplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain")    
    else:
        sns.lineplot(data=df_melt, x="RA_percentage", y="species abundance", hue="strain",orient="y")

    plt.xlabel("Percentage of maximal RA production")

    plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()


def plot_biomass_time_course(sim_results, exp_data, growth_curves = None):
    """Make biomass time-course plot for CAL11:SAL11:MAM3 triculture"""

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
    sns.lineplot(x='time', y='biomass', data=exp_data, color = "red", label="measured total biomass", linestyle ="-.", legend=False)

    if growth_curves is not None:
        sns.lineplot(x='time', y='biomass', hue="strain", data=growth_curves, linestyle ="-.", hue_order=["CAL11", "SAL11", "MAM3"], legend=False)

    plt.legend()
    plt.xlabel("time (h)")
    plt.ylabel("biomass (gDW)")


def plot_biomass_time_course_glc(sim_results, exp_data):
    """Make biomass time-course plot for CAL2:SAL9:MAM2 triculture"""

    """ sim_results should be a dataframe outputted from sim.total_biomass """

    biomass_df = sim_results.copy()
    
    # convert from cycles to time, assuming a simulation rate of 0.1 hours
    biomass_df["time"] = biomass_df["cycle"]*0.1
    biomass_df.drop(columns=["cycle"], inplace=True)

    # add total biomass as a column
    biomass_df["total_BM"] = biomass_df["CAL2"] + biomass_df["SAL9"] + biomass_df["MAM2"]

    # prepare df for plotting
    plot_df = biomass_df.melt('time', var_name='variable', value_name='biomass')

    # plot
    sns.lineplot(data=plot_df, x="time", y="biomass", hue="variable")

    plt.lineplot(x='time', y='biomass', data=exp_data, color = "red", linestyle ="-.", label="measured total biomass", legend=False)
    plt.legend()
    plt.xlabel("time (h)")
    plt.ylabel("biomass (gDW)")


def plot_biomass_time_course_coculture(sim_results):
    """Make biomass time-course plot for RAU2:RAD4 coculture"""

    """ sim_results should be a dataframe outputted from sim.total_biomass """

    biomass_df = sim_results.copy()
    
    # convert from cycles to time, assuming a simulation rate of 0.1 hours
    biomass_df["time"] = biomass_df["cycle"]*0.1
    biomass_df.drop(columns=["cycle"], inplace=True)

    # add total biomass as a column
    biomass_df["total_BM"] = biomass_df["RAU2"] + biomass_df["RAD4"]

    # prepare df for plotting
    plot_df = biomass_df.melt('time', var_name='variable', value_name='biomass')

    # plot
    sns.lineplot(data=plot_df, x="time", y="biomass", hue="variable")

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
    sns.lineplot(data=exp_data, x='time', y='subpopulation_percentage', hue='strain', hue_order=["CAL11", "SAL11", "MAM3"], linestyle ="-.", legend=False)
    
    plt.ylabel("Subpopulation fraction")
    plt.xlabel("Time (h)")


def plot_relative_abundance_time_course_coculture(sim_results):
    """Make biomass time-course plot for RAU2:RAD4 coculture"""

    """ sim_results should be a dataframe outputted from sim.total_biomass """
    
    biomass_df = sim_results.copy()

    total_BM = biomass_df["RAU2"] + biomass_df["RAD4"]
    RAU2_frac = biomass_df["RAU2"] / total_BM
    RAD4_frac = biomass_df["RAD4"] / total_BM

    df = pd.concat([total_BM, RAU2_frac, RAD4_frac], axis=1, keys=["total_BM", "RAU2", "RAD4"])
    df["time"] = df.index * 0.1
    plot_df = df.melt(id_vars="time", value_vars=["RAU2", "RAD4"], value_name="relative_abundance", var_name="strain")

    sns.lineplot(plot_df, x="time", y="relative_abundance", hue="strain", hue_order=["RAU2", "RAD4"])
    
    plt.ylabel("Subpopulation fraction")
    plt.xlabel("Time (h)")


def plot_relative_abundance_time_course_glc(sim_results, exp_data):
    """Make biomass time-course plot for CAL2:SAL9:MAM2 triculture"""

    """ sim_results should be a dataframe outputted from sim.total_biomass """
    
    biomass_df = sim_results.copy()

    total_BM = biomass_df["CAL2"] + biomass_df["SAL9"] + biomass_df["MAM2"]
    CAL2_frac = biomass_df["CAL2"] / total_BM
    SAL9_frac = biomass_df["SAL9"] / total_BM
    MAM2_frac = biomass_df["MAM2"] / total_BM

    df = pd.concat([total_BM, CAL2_frac, SAL9_frac, MAM2_frac], axis=1, keys=["total_BM", "CAL2", "SAL9", "MAM2"])
    df["time"] = df.index * 0.1
    plot_df = df.melt(id_vars="time", value_vars=["CAL2", "SAL9", "MAM2"], value_name="relative_abundance", var_name="strain")

    sns.lineplot(plot_df, x="time", y="relative_abundance", hue="strain", hue_order=["CAL2", "SAL9", "MAM2"])
    sns.lineplot(data=exp_data, x='time', y='subpopulation_percentage', hue='strain', hue_order=["CAL2", "SAL9", "MAM2"], linestyle="-.", legend=False)
    
    plt.ylabel("Subpopulation fraction")
    plt.xlabel("Time (h)")


def plot_production_time_course(sim_results, exp_data = None):

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

    if exp_data is not None:
        # create a second y-axis (ax2) and plot scatterplot on it
        ax2 = ax1.twinx()
        sns.lineplot(data=exp_data, x="time", y="mmol_per_L", hue="product", hue_order=["CA", "SAA", "RA"], ax=ax2, linestyle ="-.")
        ax2.set_ylabel("measured concentration (mmol/L)")

    plt.show()


def plot_production_flux_values(sim):
    """Plot the flux values for the strains."""

    CAL11_flux = sim.get_species_exchange_fluxes("CAL11")
    SAL11_flux = sim.get_species_exchange_fluxes("SAL11")
    MAM3_flux = sim.get_species_exchange_fluxes("MAM3")

    production_fluxes = pd.concat([CAL11_flux["EX_34dhcinm_e"], SAL11_flux["EX_saa_e"], MAM3_flux["EX_rosma_e"]], axis=1)
    production_fluxes["time"] = production_fluxes.index *0.1
    plot_df = production_fluxes.melt(id_vars="time", value_vars=["EX_34dhcinm_e", "EX_saa_e", "EX_rosma_e"], value_name="flux_value", var_name="reaction")

    sns.lineplot(data=plot_df, x="time", y="flux_value", hue="reaction")


def plot_production_flux_values_glc(sim):
    """Plot the flux values for the strains. Ids are hard-coded for the CAL2:SAL9:MAM2 triculture"""

    CAL2_flux = sim.get_species_exchange_fluxes("CAL2")
    SAL9_flux = sim.get_species_exchange_fluxes("SAL9")
    MAM2_flux = sim.get_species_exchange_fluxes("MAM2")

    production_fluxes = pd.concat([CAL2_flux["EX_34dhcinm_e"], SAL9_flux["EX_saa_e"], MAM2_flux["EX_rosma_e"]], axis=1)
    production_fluxes["time"] = production_fluxes.index *0.1
    plot_df = production_fluxes.melt(id_vars="time", value_vars=["EX_34dhcinm_e", "EX_saa_e", "EX_rosma_e"], value_name="flux_value", var_name="reaction")

    sns.lineplot(data=plot_df, x="time", y="flux_value", hue="reaction")


def plot_production_flux_values_coculture(sim):
    """Plot the flux values for the strains. Ids are hard-coded for the RAU2:RAD4 coculture"""

    RAU2_flux = sim.get_species_exchange_fluxes("RAU2")
    RAD4_flux = sim.get_species_exchange_fluxes("RAD4")

    production_fluxes = pd.concat([RAU2_flux["EX_34dhcinm_e"], RAD4_flux["EX_saa_e"], RAD4_flux["EX_rosma_e"]], axis=1)
    production_fluxes["time"] = production_fluxes.index *0.1
    plot_df = production_fluxes.melt(id_vars="time", value_vars=["EX_34dhcinm_e", "EX_saa_e", "EX_rosma_e"], value_name="flux_value", var_name="reaction")

    sns.lineplot(data=plot_df, x="time", y="flux_value", hue="reaction")


def plot_metabolites(sim):
    """Plot metabolite amounts."""
    
    sim_results = sim.get_metabolite_time_series(upper_threshold = 900.)
    sim_results.plot(x = "cycle")
    plt.ylabel("mmol")


def add_missing_vals(df):
    """add the rows that didn't compute as 0 values inplace, hard coded."""

    df = pd.concat([df, 
        pd.DataFrame([{"inoculation_ratio": "(2, 1, 1)", "glc_xyl_ratio":"(1, 4)", "total_biomass": 0, "total_RA": 0}]),
        pd.DataFrame([{"inoculation_ratio": "(3, 1, 1)", "glc_xyl_ratio":"(1, 4)", "total_biomass": 0, "total_RA": 0}])
        ], axis=0, ignore_index=True)
    

def mmol_to_mg_L(mmol, product = "RA"):
        """Helper function"""

        if product=="RA":
            MM = 360.3148
        elif product=="CA":
            MM = 180.1574
        elif product=="SAA":
            MM = 198.1727

        mg = mmol*MM
        mg_L = mg/0.1
        return mg_L


def plot_inoculum_substrate(results_df):

    df = results_df.copy()

    if results_df.shape[0] != 36:
        add_missing_vals(df) 
    
    df["RA_concentration"] = mmol_to_mg_L(df["total_RA"])

    fig_A = df[df["glc_xyl_ratio"] == "(1, 4)"]
    fig_B = df[df["glc_xyl_ratio"] == "(2, 3)"]
    fig_C = df[df["glc_xyl_ratio"] == "(3, 2)"]
    fig_D = df[df["glc_xyl_ratio"] == "(4, 1)"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 6))

    sns.barplot(data=fig_A, x="inoculation_ratio", y="RA_concentration", ax=axes[0, 0])
    axes[0, 0].set_title('xylose:glucose=4:1')

    sns.barplot(data=fig_B, x="inoculation_ratio", y="RA_concentration", ax=axes[0, 1])
    axes[0, 1].set_title('xylose:glucose=3:2')

    sns.barplot(data=fig_C, x="inoculation_ratio", y="RA_concentration", ax=axes[1, 0])
    axes[1, 0].set_title('xylose:glucose=2:3')

    sns.barplot(data=fig_D, x="inoculation_ratio", y="RA_concentration", ax=axes[1, 1])
    axes[1, 1].set_title('xylose:glucose=1:4')

    for ax in axes.flat:
        ax.set_ylabel('RA concentration (mg/L)')

    plt.tight_layout()

    plt.show()      

    
def plot_inoculum_4D(results_df):

    df = results_df.copy()

    df["RA_concentration"] = mmol_to_mg_L(df["total_RA"])

    fig, ax = plt.subplots(figsize=(10, 5))

    sns.barplot(data=df, x="inoculation_ratio", y="RA_concentration")

    ax.set_ylabel('RA concentration (mg/L)')
    ax.set_xlabel("inoculation ratio")
    plt.tight_layout()
    plt.show()   


def plot_inoculum_3D(results_df):

    df = results_df.copy()

    df["CA_concentration"] = mmol_to_mg_L(df["total_CA"], product="CA")
    df["SAA_concentration"] = mmol_to_mg_L(df["total_SAA"], product="SAA")
    df["RA_concentration"] = mmol_to_mg_L(df["total_RA"], product="RA")

    # drop columns
    drop_cols = ['total_biomass', 'total_RA', 'total_CA','total_SAA']
    df.drop(columns=drop_cols, inplace=True)
    
    # rename coumns
    df.rename(columns={"CA_concentration": "CA", "SAA_concentration": "SAA", "RA_concentration": "RA"}, inplace=True)

    # melt df on CA, SAA, RA
    plot_df = df.melt(id_vars="inoculation_ratio", value_vars=["CA","SAA", "RA"], value_name="concentration", var_name="product")

    fig, ax = plt.subplots(figsize=(10, 5))

    sns.barplot(data=plot_df, x="inoculation_ratio", y="concentration", hue="product")

    ax.set_ylabel('product concentration (mg/L)')
    ax.set_xlabel("inoculation ratio")
    plt.tight_layout()
    plt.show()  


def get_best_performing_combo(df):
    df_with_RA = df.copy()
    df_with_RA["RA_concentration"] = mmol_to_mg_L(df_with_RA["total_RA"])
    print("The highest RA production with the combination:\n", df_with_RA.loc[df_with_RA['RA_concentration'].idxmax()])


def save_fig(filepath):
    """Save figure as vector image"""
    #TODO: implement
    pass