"""
Running dFBA using COMETSpy.

TODO: create a better soluton for the MM kinetic paramerers
"""

import cometspy as c
from cobra.io import read_sbml_model
import matplotlib.pyplot as plt
import pandas as pd
import warnings


def simulate_xyl_glc_triculture(cal11, sal11, mam3, initial_pop_ratio: tuple[int] =(2, 1, 1), 
                                adjust_atp_requirements: bool=False, RA_lb: float = 0.00104, 
                                glc_xyl_mmol: tuple[float] = (1.67, 1.33), initial_pop: float = None) -> c.comets:

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
    if not initial_pop:
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
    #NOTE: changing up these params to see what happens - if this does not have an effect, try changing after params are set? (move params up)

    glc_vmax_adjustment = 3
    xyl_vmax_adjustment = 0.55
    km_xyl_adj = 75
    km_glc_adj = 5500

    O2Vmax = 15
    O2Km = 0.024 #mmol/L
    GlcVmax = 10.5 * glc_vmax_adjustment
    GlcKm = 0.000015  * km_glc_adj # mmol/ml #0.0027 g/L
    XylVmax = 6 * xyl_vmax_adjustment
    XylKm = 0.00011 * km_xyl_adj # mmol/ml #0.0165 g/L

    cal11_c.change_vmax("EX_o2_e", O2Vmax)
    cal11_c.change_km("EX_o2_e", O2Km)
    cal11_c.change_vmax("EX_xyl__D_e", XylVmax)
    cal11_c.change_km("EX_xyl__D_e", XylKm)

    sal11_c.change_vmax("EX_o2_e", O2Vmax)
    sal11_c.change_km("EX_o2_e", O2Km)
    sal11_c.change_vmax("EX_glc__D_e", GlcVmax)
    sal11_c.change_km("EX_glc__D_e", GlcKm)

    mam3_c.change_vmax("EX_o2_e", O2Vmax)
    mam3_c.change_km("EX_o2_e", O2Km)
    mam3_c.change_vmax("EX_xyl__D_e", XylVmax)
    mam3_c.change_km("EX_xyl__D_e", XylKm)

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


def simulate_glc_triculture(cal2, sal9, mam2, initial_pop_ratio: tuple[int] =(2, 3, 1), 
                                adjust_atp_requirements: bool=False, RA_lb: float = 0.00104, 
                                glc_mmol: float = 2.78, initial_pop: float = 1.e-3, glc_vmax_adjustment = 1, km_glc_adj = 1) -> c.comets:

    if adjust_atp_requirements:
        # in order to equal the playing field between BL21 and K12 derived models:

        # adjust the biomass reaction
        sal9.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"atp_c":-75.55223}, combine=False)
        sal9.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"h_c":75.377230}, combine=False)
        sal9.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"adp_c":75.377230}, combine=False)
        sal9.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"pi_c":75.373230}, combine=False)
        sal9.reactions.get_by_id("BIOMASS_Ec_iHK1487_core").add_metabolites({"h2o_c":-70.028756}, combine=False)

        # adjust the ATP maintanance requrenment
        sal9.reactions.ATPM.lower_bound = 6.86 
    
    # make comets models
    cal2_c = c.model(cal2)
    sal9_c = c.model(sal9)
    mam2_c = c.model(mam2)

    cal2_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[0]]
    sal9_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[1]]
    mam2_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[2]]

    # open the exhange reactions
    cal2_c.open_exchanges()
    sal9_c.open_exchanges()
    mam2_c.open_exchanges()

    # set to run pFBA
    cal2_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    sal9_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    mam2_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"

    # set MM kinetic parameters for glucose, oxygen, and xylose uptake reactions

    # glc_vmax_adjustment = 4
    # km_glc_adj = 10000

    O2Vmax = 15
    O2Km = 0.024 #mmol/L
    GlcVmax = 10.5 * glc_vmax_adjustment
    GlcKm = 0.000015  * km_glc_adj # mmol/ml #0.0027 g/L


    cal2_c.change_vmax("EX_o2_e", O2Vmax)
    cal2_c.change_km("EX_o2_e", O2Km)
    cal2_c.change_vmax("EX_glc__D_e", GlcVmax)
    cal2_c.change_km("EX_glc__D_e", GlcKm)

    sal9_c.change_vmax("EX_o2_e", O2Vmax)
    sal9_c.change_km("EX_o2_e", O2Km)
    sal9_c.change_vmax("EX_glc__D_e", GlcVmax)
    sal9_c.change_km("EX_glc__D_e", GlcKm)

    mam2_c.change_vmax("EX_o2_e", O2Vmax)
    mam2_c.change_km("EX_o2_e", O2Km)
    mam2_c.change_vmax("EX_glc__D_e", GlcVmax)
    mam2_c.change_km("EX_glc__D_e", GlcKm)

    # create a 1x1 layout

    layout = c.layout([cal2_c, sal9_c, mam2_c])

    # set metabolite availability

    unlimited_mets = ['ca2_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'h2o_e', 'k_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e','o2_e', 'pi_e', 'sel_e', 'slnt_e', 'so4_e', 'tungs_e', 'zn2_e']

    for met in unlimited_mets:
        layout.set_specific_metabolite(met, 1000.)
        
    layout.set_specific_metabolite("glc__D_e",  glc_mmol)
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
    params.set_param("maxCycles", 500) # 0.1 hours x 500 cycles for 50 hours simulated

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


def simulate_coculture(rau2, rad4, initial_pop_ratio: tuple[int] =(3, 1), RA_lb: float = 0.00104, 
                       glc_mmol: float = 2.78, initial_pop: float = 1.e-3, glc_vmax_adjustment = 1, km_glc_adj = 1) -> c.comets:

    
    # make comets models
    rau2_c = c.model(rau2)
    rad4_c = c.model(rad4)

    rau2_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[0]]
    rad4_c.initial_pop = [0,0,initial_pop*initial_pop_ratio[1]]

    # open the exhange reactions
    rau2_c.open_exchanges()
    rad4_c.open_exchanges()

    # set to run pFBA
    rau2_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"
    rad4_c.obj_style="MAX_OBJECTIVE_MIN_TOTAL"

    # set MM kinetic parameters for glucose, oxygen, and xylose uptake reactions

    O2Vmax = 15
    O2Km = 0.024 #mmol/L
    GlcVmax = 10.5 * glc_vmax_adjustment
    GlcKm = 0.000015  * km_glc_adj # mmol/ml #0.0027 g/L

    rau2_c.change_vmax("EX_o2_e", O2Vmax)
    rau2_c.change_km("EX_o2_e", O2Km)
    rau2_c.change_vmax("EX_glc__D_e", GlcVmax)
    rau2_c.change_km("EX_glc__D_e", GlcKm)

    rad4_c.change_vmax("EX_o2_e", O2Vmax)
    rad4_c.change_km("EX_o2_e", O2Km)
    rad4_c.change_vmax("EX_glc__D_e", GlcVmax)
    rad4_c.change_km("EX_glc__D_e", GlcKm)

    # create a 1x1 layout

    layout = c.layout([rau2_c, rad4_c])

    # set metabolite availability

    unlimited_mets = ['ca2_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e','h_e', 'h2o_e', 'k_e', 'mg2_e', 
                    'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e','o2_e', 'pi_e', 'sel_e', 'slnt_e', 'so4_e', 'tungs_e', 'zn2_e']

    for met in unlimited_mets:
        layout.set_specific_metabolite(met, 1000.)
        
    layout.set_specific_metabolite("glc__D_e",  glc_mmol)
    layout.set_specific_metabolite("phe__L_e", 0.05)
    layout.set_specific_metabolite("tyr__L_e", 0.05)

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
    params.set_param("maxCycles", 480) # 0.1 hours x 500 cycles for 50 hours simulated

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


def run_dfba(initial_pop_ratio: tuple[int] =(2, 1, 1), adjust_atp_requirements: bool=False, 
             RA_lb: float = 0.00104, production_lbs: tuple[float] = (0.00104, 0.00104, 0.00104), 
             glc_xyl_mmol: tuple[float] = (1.67, 1.33)) -> c.comets:
    
    """production_lbs should be in order cal, sal, mam"""

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
    mam3.reactions.RAt.bounds = (production_lbs[2], 1000.0)
    sal11.reactions.SAAt.bounds = (production_lbs[1], 1000.0)
    cal11.reactions.get_by_id("34DHCINMt").bounds = (production_lbs[0], 1000.0)
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