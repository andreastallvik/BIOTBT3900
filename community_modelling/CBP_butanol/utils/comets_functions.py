"""Run dFBA simulations using COMETS."""

import cometspy as c
import warnings
import seaborn as sns


UNLIMITED_METABOLITES = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'k_e', 'h2o_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'pi_e', 'so4_e', 'zn2_e']

SPACE_WIDTH = 3.684


def single_strain(model, medium: dict = {}, initial_pop: float = 1.e-3, sim_time: float = 140):

    # make a comets model
    comets_model = c.model(model)

    # set initial population
    comets_model.initial_pop = [0, 0, initial_pop]

    # open all exhange reactions
    comets_model.open_exchanges()

    # use pFBA when solving
    comets_model.obj_style="MAX_OBJECTIVE_MIN_TOTAL"

    # TODO: Set MM kinetics for uptake reactions

    # create a 1x1 layout
    layout = c.layout([comets_model])

    # set metabolite availability
    for met in UNLIMITED_METABOLITES:
        layout.set_specific_metabolite(met, 1000.)
    
    for met, mmol in medium.items():
        layout.set_specific_metabolite(met, mmol)

    # create params object
    params = c.params()

    # set grid size specifications
    params.set_param("spaceWidth", SPACE_WIDTH)

    # set simulation parameters
    time_step = 0.1 # hours
    params.set_param("timeStep", time_step)
    params.set_param("maxSpaceBiomass", 10.) # max gDW in simulation TODO: set this to a sensible number and not an arbitrarily large one
    params.set_param("maxCycles", int(sim_time / time_step))

    # set logging parameters
    params.set_param("writeFluxLog", True)
    params.set_param("writeMediaLog", True)
    params.set_param("FluxLogRate", 1)
    params.set_param("MediaLogRate", 1)

    # create simultion object
    sim = c.comets(layout, params)

    with warnings.catch_warnings(): #to avoid getting many futurewarning messages from a cometspy function
        warnings.simplefilter(action='ignore', category=FutureWarning)
        sim.run()

    return sim


def plot_metabolites(sim, metabolites, time_step = 0.1):
    """Plot specific metabolites from comets simulation results."""

    # get metabolite time-course data
    metabolites_time_series = sim.get_metabolite_time_series()

    # prepare dataframe for plotting

    # retreieve columns with desired metabolites
    time = metabolites_time_series["cycle"] * time_step
    present_metabolites = [m for m in metabolites if m in metabolites_time_series.columns]
    df = metabolites_time_series[present_metabolites].copy()
    missing_metabolites = list(set(metabolites) - set(present_metabolites))
    for m in missing_metabolites:
        df[m] = 0

    # convert mmol to g/L
    df = df.apply(lambda x: mmol_to_g_per_L(x.name, x))
    
    # add time column
    df["time"] = time

    plot_df = df.melt(id_vars="time", value_name="g/L")
    
    sns.lineplot(data=plot_df, x="time", y="g/L", hue="metabolite")


def plot_biomass(sim, time_step=0.1):
    """Plot biomass from comets simulation results."""

    biomass_time_series = sim.total_biomass.copy()
    time = biomass_time_series["cycle"] * time_step
    biomass_time_series["time"] = time
    biomass_time_series.drop(columns = ["cycle"], inplace=True)
    plot_df = biomass_time_series.melt(id_vars="time", var_name="strain", value_name="biomass")
    sns.lineplot(data=plot_df, x="time", y="biomass", hue="strain")


def mmol_to_g_per_L(met_name, met_mmol, volume = 0.05):
    """Convert mmol to g/L.""" 
    
    # dictionary of molar mass for each metabolite
    MM = {"xyl__D_e": 150.13, 
          "etoh_e": 46.068, 
          "but_e": 88.11, 
          "btoh_e": 74.12, 
          "ac_e": 59.044, 
          "acetone_e": 58.08}
    
    # divide by 1000 (-> mol), multiply by molar mass (-> g) and divide by volume (-> g/L)
    return (met_mmol / 1000) * MM[met_name] / volume