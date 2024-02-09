"""Run dFBA simulations using COMETS."""

import cometspy as c
import warnings
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import traceback


UNLIMITED_METABOLITES = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'k_e', 'h2o_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'pi_e', 'so4_e', 'zn2_e']

SPACE_WIDTH = 3.684


def single_strain(model, medium: dict = {}, initial_pop: float = 1.e-3, sim_time: float = 140, km_dict: dict = {}, vmax_dict: dict = {}):
    """Run a comets simulation for a single strain

    Args:
        model: cobrapy model
        medium (dict, optional): metabolite: molar amount to be included in sim. Defaults to {}. If empty, only std. unlimited metabolites will be included.
        initial_pop (float, optional): initial biomass of organism. Defaults to 1.e-3.
        sim_time (float, optional): number of hours in simulation. Defaults to 140.
        km_dict (dict, optional): dict of rx:value for km values to set. Defaults to {}.
        vmax_dict (dict, optional): dict of rx:value for vmax values to set. Defaults to {}.

    Returns:
        sim: comets simulation object
    """

    # make a comets model
    comets_model = c.model(model)

    # set initial population
    comets_model.initial_pop = [0, 0, initial_pop]

    # open all exhange reactions
    comets_model.open_exchanges()

    # use pFBA when solving
    comets_model.obj_style="MAX_OBJECTIVE_MIN_TOTAL"

    # set MM kinetic parameters
    set_kinetic_params(comets_model, vmax_dict, km_dict)

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

    try:
        with warnings.catch_warnings(): #to avoid getting many futurewarning messages from a cometspy function
            warnings.simplefilter(action='ignore', category=FutureWarning)
            sim.run()

    except Exception as e:
        with open('run_output.log', 'w') as f:
            f.write(sim.run_output)
            f.write('\n\nException traceback:\n')
            f.write(traceback.format_exc())
        raise  # re-raise the exception after logging

    return sim


def mult_strain(models: list, medium: dict = {}, initial_pop: float = 1.e-3, sim_time: float = 140, specific_initial_pop: dict = {}, kinetic_params: dict = {}):
    """Run a simulation for multiple strains.

    Args:
        models (list): list of cobrapy models
        medium (dict, optional): Dict of metabolite:molar amount, if empty only unlimited metabolites are included. Defaults to {}.
        initial_pop (float, optional): initial biomass for all community members, ignored if specific_initial_pop is non-empty. Defaults to 1.e-3.
        sim_time (float, optional): number of hours of simulation time. Defaults to 140.
        specific_initial_pop (dict, optional): dictionary of cobrapy model.id:initial biomass. If not empty this is used instead of initial_pop. Defaults to {}.
        kinetic_params (dict, optional): dict of model_id:{"vmax":{rx:val}, "km":{rx:val}} for setting kinetic parameters for each model. Defaults to {}.

    Returns:
        sim: comets simulation object
    """

    comets_models = {model.id:None for model in models}

    # for each model

    for model in models:

        # make a comets model
        comets_models[model.id] = c.model(model)

        # set initial population
        if specific_initial_pop:
            initial_pop = specific_initial_pop[model.id]
            
        comets_models[model.id].initial_pop = [0, 0, initial_pop]
        
        # open all exhange reactions
        comets_models[model.id].open_exchanges()

        # use pFBA when solving
        comets_models[model.id].obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    

    # set kinetic parameters
    for model_id, params in kinetic_params.items():
        vmax_dict = params.get("vmax", {})
        km_dict = params.get("km", {})
        set_kinetic_params(comets_models[model_id], vmax_dict, km_dict)

    # create a 1x1 layout
    layout = c.layout(list(comets_models.values()))

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


def sequental_com(m5, nj4, init_medium: dict = {}, initial_pop_m5: float = 1.e-3, initial_pop_nj4: float = 1.e-3, 
                  total_sim_time: float = 140, inoc_time: float = 60, kinetic_params: dict = {}):
    """Run 2 sequential COMETS simulations for the CBP butanol community.

    Args:
        m5 (cobrapy model): GEM for m5
        nj4 (cobrapy model): GEM for nj4
        init_medium (dict, optional): Culture medium. If empty, only UNLIMITED METABOLITES are added. Therefore carbon source must be added. Defaults to {}.
        initial_pop_m5 (float, optional): Initail biomass for m5 strain. Defaults to 1.e-3.
        initial_pop_nj4 (float, optional): Initial biomass for nj4 strain. Defaults to 1.e-3.
        total_sim_time (float, optional): Total hours to simulate for. Defaults to 140.
        inoc_time (float, optional): Hour in which nj4 is added to the community. Defaults to 60.
        kinetic_params (dict, optional): dict of model_id:{"vmax":{rx:val}, "km":{rx:val}} for setting kinetic parameters for each model. Defaults to {}.

    Returns:
        tuple(c.sim, c.sim): Simulation object for the first sim (only m5) and the second sim (m5 + nj4).
    """
    

    #TODO: add in some failsafes for if the simulation fails
    
    second_sim_time = total_sim_time - inoc_time

    # get the kinetic parameters for the first sim if any
    m5_kinetic_params = kinetic_params.get(m5.id, {})
    m5_vmax = m5_kinetic_params.get("vmax", {})
    m5_km = m5_kinetic_params.get("km", {})

    # run a single-strain simulation for the first strain
    first_sim = single_strain(m5, medium=init_medium, initial_pop=initial_pop_m5, sim_time=inoc_time, km_dict=m5_km, vmax_dict=m5_vmax)

    # retrieve information from the first simulation 
    
    biomass_m5 = first_sim.total_biomass["M5"].iloc[-1]
    
    # final metabolite amounts in the medium
    metabolites = first_sim.get_metabolite_time_series().iloc[-1, 2:]
    new_medium = {met:mol for met,mol in metabolites.items() if mol > 0.0}

    # run a mult-strain simulation for the second strain
    second_sim = mult_strain([m5, nj4], medium=new_medium, sim_time=second_sim_time, specific_initial_pop={"NJ4":initial_pop_nj4, "M5": biomass_m5}, kinetic_params=kinetic_params)

    return first_sim, second_sim


def set_kinetic_params(model: c.model, vmax_dict: dict = {}, km_dict: dict = {}):
    """Set kinetic parameters for a comets model inplace.

    Args:
        model (c.model): comets model object
        vmax_dict (dict, optional): dict of rx:value for vmax values to set. Defaults to {}.
        km_dict (dict, optional): dict of rx:value for km values to set. Defaults to {}.
    """
    
    for rx, vmax in vmax_dict.items():
        model.change_vmax(rx, vmax)

    for rx, km in km_dict.items():
        model.change_km(rx, km)


def collapse_sequential_sim(sim_1, sim_2):
    """Helper function to combine the results from 2 sequential simulation objects generated by sequential_com() into dataframes for analysis.

    Args:
        sim_1 (c.sim): first simulation object
        sim_2 (c.sim): second simulation object

    Returns:
        df: dataframe of biomass, metabolite and a dictionary of flux time-series data for each strain.
    """

    # biomass
    bm_1 = sim_1.total_biomass.copy()
    bm_2 = sim_2.total_biomass.copy()

    cycle_diff = bm_1["cycle"].iloc[bm_1.shape[0]-1] + 1

    bm_2["cycle"] = bm_2["cycle"] + cycle_diff
    bm = pd.concat([bm_1, bm_2]).fillna(value=0)

    # metabolies
    met_1 = sim_1.get_metabolite_time_series()
    met_2 = sim_2.get_metabolite_time_series()

    met_2["cycle"] = met_2["cycle"] + cycle_diff
    met = pd.concat([met_1, met_2]).fillna(value=0)

    # fluxes
    # M5 fluxes
    flux_1_a = sim_1.fluxes_by_species["M5"].copy()
    flux_2_a = sim_2.fluxes_by_species["M5"].copy()
    flux_2_a["cycle"] = flux_2_a["cycle"] + cycle_diff
    flux_1 = pd.concat([flux_1_a, flux_2_a]).fillna(value=0)

    #NJ4 fluxes
    flux_2 = sim_2.fluxes_by_species["NJ4"].copy()
    flux_2["cycle"] = flux_2["cycle"] + cycle_diff
    # adding in 0s for the cycles before NJ4 is added
    zeros_df = pd.DataFrame(0, index=range(cycle_diff), columns=flux_2.columns)
    zeros_df["cycle"] = np.arange(1,cycle_diff+1)
    flux_2 = pd.concat([zeros_df, flux_2], ignore_index=True)

    fluxes = {"M5": flux_1, "NJ4": flux_2}

    return bm, met, fluxes


def plot_metabolites(sim = None, metabolites = None, time_step = 0.1, metabolites_time_series = None, inoc_time = None):
    """Plot specific metabolites from comets simulation results."""

    if sim is None and metabolites_time_series is None:
        raise ValueError("Either a comets simulation object or a dataframe of metabolite time-series data must be provided.")

    if metabolites is None:
        raise ValueError("A list of metabolites to plot must be provided.")

    # get metabolite time-course data
    if metabolites_time_series is None:
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

    if inoc_time is not None:
        plt.axvline(x=inoc_time, color='k', linestyle='--')


def plot_biomass(sim = None, time_step=0.1, total_biomass = None, inoc_time = None):
    """Plot biomass from comets simulation results."""

    if sim is None and total_biomass is None:
        raise ValueError("Either a comets simulation object or a dataframe of biomass time-series data must be provided.")

    if total_biomass is None:
        biomass_time_series = sim.total_biomass.copy()
    else:
        biomass_time_series = total_biomass.copy()

    time = biomass_time_series["cycle"] * time_step
    biomass_time_series["time"] = time
    biomass_time_series.drop(columns = ["cycle"], inplace=True)
    plot_df = biomass_time_series.melt(id_vars="time", var_name="strain", value_name="biomass")
    sns.lineplot(data=plot_df, x="time", y="biomass", hue="strain")
    if inoc_time is not None:
        plt.axvline(x=inoc_time, color='k', linestyle='--')
        

def mmol_to_g_per_L(met_name, met_mmol, volume = 0.05):
    """Convert mmol to g/L.""" 
    
    # dictionary of molar mass for each metabolite
    MM = {"xyl__D_e": 150.13, 
          "etoh_e": 46.068, 
          "but_e": 88.11, 
          "btoh_e": 74.12, 
          "ac_e": 59.044, 
          "acetone_e": 58.08,
          "xylan4_e": 600.52,
          "xylan8_e": 1201.04}
    
    # divide by 1000 (-> mol), multiply by molar mass (-> g) and divide by volume (-> g/L)
    return (met_mmol / 1000) * MM[met_name] / volume


def plot_reaction_flux(sim = None, reactions: list = None, strain: str = None, time_step=0.1, fluxes = None, inoc_time=None):
    """Plot reaction fluxes for a single strain from a comets simulation."""
    
    if sim is None and fluxes is None:
        raise ValueError("Either a comets simulation object or a dataframe of fluxes time-series data must be provided.")

    if reactions is None:
        raise ValueError("A list of reactions to plot must be provided.")
    
    if fluxes is None:
        fluxes = sim.fluxes_by_species[strain]

    time = fluxes["cycle"] * time_step

    present_reaction_fluxes = [rx for rx in reactions if rx in fluxes.columns]
    df = fluxes[present_reaction_fluxes].copy()

    missing_reactions = list(set(reactions) - set(present_reaction_fluxes))
    for rx in missing_reactions:
            df[rx] = 0
    df["time"] = time

    plot_df = df.melt(id_vars="time", value_name="flux", var_name="reaction")
    
    sns.lineplot(data=plot_df, x="time", y="flux", hue="reaction") 

    if inoc_time is not None:
        plt.axvline(x=inoc_time, color='k', linestyle='--')
        