"""Run dFBA simulations using COMETS."""

import cometspy as c
import warnings


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
    max_cycles = int(sim_time / time_step) #TODO: write in some rounding thing here to ensure no bugs
    print(max_cycles)
    params.set_param("timeStep", time_step)
    params.set_param("maxSpaceBiomass", 10.) # max gDW in simulation TODO: set this to a sensible number and not an arbitrarily large one
    params.set_param("maxCycles", max_cycles)

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

