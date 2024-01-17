"""Run dFBA simulations using COMETS."""

import cometspy as c
import warnings
import seaborn as sns


UNLIMITED_METABOLITES = ['ca2_e', 'cl_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'k_e', 'h2o_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'pi_e', 'so4_e', 'zn2_e']

SPACE_WIDTH = 3.684


def single_strain(model, medium: dict = {}, initial_pop: float = 1.e-3, sim_time: float = 140):
    """Run a comets simulation for a single strain

    Args:
        model: cobrapy model
        medium (dict, optional): metabolite: molar amount to be included in sim. Defaults to {}. If empty, only std. unlimited metabolites will be included.
        initial_pop (float, optional): initial biomass of organism. Defaults to 1.e-3.
        sim_time (float, optional): number of hours in simulation. Defaults to 140.

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


def mult_strain(models: list, medium: dict = {}, initial_pop: float = 1.e-3, sim_time: float = 140, specific_initial_pop: dict = {}):
    """Run a simulation for multiple strains.

    Args:
        models (list): list of cobrapy models
        medium (dict, optional): Dict of metabolite:molar amount, if empty only unlimited metabolites are included. Defaults to {}.
        initial_pop (float, optional): initial biomass for all community members, ignored if specific_initial_pop is non-empty. Defaults to 1.e-3.
        sim_time (float, optional): number of hours of simulation time. Defaults to 140.
        specific_initial_pop (dict, optional): dictionary of cobrapy model.id:initial biomass. If not empty this is used instead of initial_pop. Defaults to {}.

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
    

    # TODO: Set MM kinetics for uptake reactions

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
                  total_sim_time: float = 140, inoc_time: float = 60):
    

    #TODO: check this works / add in some failsafes for if the simulation fails
    
    second_sim_time = total_sim_time - inoc_time

    # run a single-strain simulation for the first strain
    first_sim = single_strain(m5, medium=init_medium, initial_pop=initial_pop_m5, sim_time=inoc_time)

    # retrieve information from the first simulation 
    
    biomass_m5 = first_sim.total_biomass["M5"].iloc[-1]
    
    # final metabolite amounts in the medium
    metabolites = first_sim.get_metabolite_time_series().iloc[-1, 2:]
    new_medium = {met:mol for met,mol in metabolites.items() if mol > 0.0}

    # run a mult-strain simulation for the second strain
    second_sim = mult_strain([m5, nj4], medium=new_medium, sim_time=second_sim_time, specific_initial_pop={"NJ4":initial_pop_nj4, "M5": biomass_m5})

    # TODO: combine the reults from the first and second sim in sime meaningful way
    return first_sim, second_sim


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
          "acetone_e": 58.08,
          "xylan4_e": 600.52,
          "xylan8_e": 1201.04}
    
    # divide by 1000 (-> mol), multiply by molar mass (-> g) and divide by volume (-> g/L)
    return (met_mmol / 1000) * MM[met_name] / volume


def plot_reaction_flux(sim, reactions: list, strain: str, time_step=0.1):
    """Plot reaction fluxes for a single strain from a comets simulation."""
    
    fluxes = sim.get_species_exchange_fluxes(strain)
    time = fluxes["cycle"] * time_step

    present_reaction_fluxes = [rx for rx in reactions if rx in fluxes.columns]
    df = fluxes[present_reaction_fluxes].copy()

    missing_reactions = list(set(reactions) - set(present_reaction_fluxes))
    for rx in missing_reactions:
            df[rx] = 0
    df["time"] = time

    plot_df = df.melt(id_vars="time", value_name="flux", var_name="reaction")
    
    sns.lineplot(data=plot_df, x="time", y="flux", hue="reaction") 