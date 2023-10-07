"""
Running dFBA using COMETSpy.

TODO: add adjustments possible through the arguments:
- different levels of glucose / xylose

"""

import cometspy as c
from cobra.io import read_sbml_model
from cobra.manipulation import knock_out_model_genes
import matplotlib.pyplot as plt
import pandas as pd
import warnings


def run_dfba(initial_pop_ratio: tuple[int] =(2, 1, 1), adjust_atp_requirements: bool=False, 
             RA_lb: float = 0.00104, glc_xyl_mmol: tuple[float] = (1.66, 1.33)) -> c.comets:

    # load models
    cal11 = read_sbml_model("../GEMs/CAL2.xml")
    sal11 = read_sbml_model("../GEMs/SAL9.xml")
    mam3 = read_sbml_model("../GEMs/MAM2.xml")

    # perform knockouts of glucose / xylose pathways
    cal11.reactions.GLCtex_copy1.bounds = (0.0, 0.0)
    cal11.reactions.GLCtex_copy2.bounds = (0.0, 0.0)
    mam3.reactions.GLCtex_copy1.bounds = (0.0, 0.0)
    mam3.reactions.GLCtex_copy2.bounds = (0.0, 0.0)
    sal11.reactions.XYLtex.bounds = (0.0, 0.0)

    # force RA/CA/SAA production at 80% of SteadyCom optimum
    mam3.reactions.RAt.bounds = (RA_lb, 1000.0)
    sal11.reactions.SAAt.bounds = (RA_lb, 1000.0)
    cal11.reactions.get_by_id("34DHCINMt").bounds = (RA_lb, 1000.0)
    # updated, using individual organisms gDW and not communty
    # mam3.reactions.RAt.bounds = (0.89, 1000.0)
    # sal11.reactions.SAAt.bounds = (0.001, 1000.0)
    # cal11.reactions.get_by_id("34DHCINMt").bounds = (0.15, 1000.0)

    # update ids
    cal11.id = "CAL11"
    sal11.id = "SAL11"
    mam3.id = "MAM3"

    if adjust_atp_requirements:
        # in order to equal the playing field between BL21 and K12 derived models:

        # adjust the biomass reaction
        sal11.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"atp_c":-75.55223}, combine=False)
        sal11.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"h_c":75.377230}, combine=False)
        sal11.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"adp_c":75.377230}, combine=False)
        sal11.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"pi_c":75.373230}, combine=False)
        sal11.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"h2o_c":-70.028756}, combine=False)

        # adjust the ATP maintanance requrenment
        sal11.reactions.ATPM.lower_bound = 6.86 
    
    # make comets models
    cal11_c = c.model(cal11)
    sal11_c = c.model(sal11)
    mam3_c = c.model(mam3)

    # set initial population
    initial_pop = 1.e-3

    cal11_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[0]]
    sal11_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[1]]
    mam3_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[2]]

    # open the exhange reactions
    cal11_c.open_exchanges()
    sal11_c.open_exchanges()
    mam3_c.open_exchanges()

    # set to run pFBA
    cal11_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    sal11_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    mam3_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"

    # set MM kinetic parameters for glucose, oxygen, and xylose uptake reactions

    cal11_c.change_vmax("EX_o2_e", 15)
    cal11_c.change_km("EX_o2_e", 0.024)
    cal11_c.change_vmax("EX_glc__D_e", 10.5)
    cal11_c.change_km("EX_glc__D_e", 0.0027)

    sal11_c.change_vmax("EX_o2_e", 15)
    sal11_c.change_km("EX_o2_e", 0.024)
    sal11_c.change_vmax("EX_xyl__D_e", 6)
    sal11_c.change_km("EX_xyl__D_e", 0.0165)

    mam3_c.change_vmax("EX_o2_e", 15)
    mam3_c.change_km("EX_o2_e", 0.024)
    mam3_c.change_vmax("EX_glc__D_e", 10.5)
    mam3_c.change_km("EX_glc__D_e", 0.0027)

    # create a 1x1 layout

    layout = c.layout([cal11_c, sal11_c, mam3_c])

    # set metabolite availability

    unlimited_mets = ['ca2_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'h2o_e', 'k_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e','o2_e', 'pi_e', 'sel_e', 'slnt_e', 'so4_e', 'tungs_e', 'zn2_e']

    for met in unlimited_mets:
        layout.set_specific_metabolite(met, 1000.)
        
    layout.set_specific_metabolite("glc__D_e",  glc_xyl_mmol[0])
    layout.set_specific_metabolite("xyl__D_e", glc_xyl_mmol[1])
    layout.set_specific_metabolite("phe__L_e", 0.05)

    #add the products so that they can be monitored
    layout.set_specific_metabolite("saa_e", 0)
    layout.set_specific_metabolite("34dhcinm_e", 0)
    layout.set_specific_metabolite("rosma_e", 0)

    # set refresh of metabolites (chemostat) #TODO: just remove this code, it doesn not do anything
    dilution_rate = 0 # / hr
    for met in unlimited_mets:
        layout.set_specific_refresh(met, 1000. * dilution_rate) # 100 mmol / hour
    layout.set_specific_refresh("glc__D_e", 1. * dilution_rate) # 0.1 mmol / hour
    layout.set_specific_refresh("xyl__D_e", 1. * dilution_rate) # 0.1 mmol / hour

    # create params object

    params = c.params()

    # set cell death and metabolite dilution
    params.set_param("deathRate", dilution_rate)
    params.set_param("metaboliteDilutionRate", dilution_rate)

    # set grid size specifications
    params.set_param("spaceWidth", 4.65)

    # set simulation parameters
    params.set_param("timeStep", 0.1) # hours
    params.set_param("maxSpaceBiomass", 10.)
    params.set_param("maxCycles", 600) # 0.1 hours x 600 cycles for 60 hours simulated

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



def run_dfba_for_monoculture(RA_lb = 0.00104) -> c.comets:

    # load models
    mra = read_sbml_model("../GEMs/MRA.xml")

    # force RA production at 80% of SteadyCom optimum
    mra.reactions.RAt.bounds = (RA_lb, 1000.0)
   
    # update ids
    mra.id = "MRA"

    
    # make comets model
    mra_c = c.model(mra)

    # set initial population
    initial_pop = 1.e-3

    mra_c.initial_pop = [0,0,initial_pop]

    # open the exhange reaction
    mra_c.open_exchanges()

    # set to run pFBA
    mra_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    
    # set MM kinetic parameters for glucose and oxygen
    mra_c.change_vmax("EX_o2_e", 15)
    mra_c.change_km("EX_o2_e", 0.024)
    mra_c.change_vmax("EX_glc__D_e", 10.5)
    mra_c.change_km("EX_glc__D_e", 0.0027)

    # create a 1x1 layout
    layout = c.layout([mra_c])

    # set metabolite availability
    unlimited_mets = ['ca2_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'h2o_e', 'k_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e','o2_e', 'pi_e', 'sel_e', 'slnt_e', 'so4_e', 'tungs_e', 'zn2_e']

    for met in unlimited_mets:
        layout.set_specific_metabolite(met, 1000.)
        
    layout.set_specific_metabolite("glc__D_e",  2.78)
    layout.set_specific_metabolite("phe__L_e", 0.05)

    #add the products so that they can be monitored
    layout.set_specific_metabolite("saa_e", 0)
    layout.set_specific_metabolite("34dhcinm_e", 0)
    layout.set_specific_metabolite("rosma_e", 0)

    # set refresh of metabolites (chemostat) #TODO: just remove this code, it doesn not do anything
    dilution_rate = 0 # / hr
    for met in unlimited_mets:
        layout.set_specific_refresh(met, 1000. * dilution_rate) # 100 mmol / hour
    layout.set_specific_refresh("glc__D_e", 1. * dilution_rate) # 0.1 mmol / hour
    layout.set_specific_refresh("xyl__D_e", 1. * dilution_rate) # 0.1 mmol / hour

    # create params object

    params = c.params()

    # set cell death and metabolite dilution
    params.set_param("deathRate", dilution_rate)
    params.set_param("metaboliteDilutionRate", dilution_rate)

    # set grid size specifications
    params.set_param("spaceWidth", 4.65)

    # set simulation parameters
    params.set_param("timeStep", 0.1) # hours
    params.set_param("maxSpaceBiomass", 10.)
    params.set_param("maxCycles", 600) # 0.1 hours x 600 cycles for 60 hours simulated

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