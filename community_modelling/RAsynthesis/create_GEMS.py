"""
Scripts to modify k12 and bl21 GEMs to equal the strains in the RA study.
"""

from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.manipulation import knock_out_model_genes
from GEM_eval import energy_generating_cycle_test
from config import ROOT_DIR
from multiprocessing import freeze_support

def no_energy_gen_cycles(model):
    """Tests if there is an energy generating cycle in the model, returns True if none is found."""
    sol = energy_generating_cycle_test(model)
    return (sol.objective_value == 0.0)
    

def production_possible(model, product_reaction):
    """Tests that the model can produce some of the required product, returns True if prodiuct is produced at 95% FVA."""
    fva_res = flux_variability_analysis(model, fraction_of_optimum=0.95)
    return fva_res["maximum"][product_reaction] > 0.0


def mam2():

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "iML1515.xml"))

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    MAM2 = k12.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = 'c'

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = 'e'

    #caffeoyl CoA
    caffcoa_c = uni.metabolites.get_by_id("caffcoa_c")
    caffcoa_c.compartment = 'c'

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #rosmarinic acid
    rosma_c = Metabolite(
        'rosma_c',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='c')

    #external rosmarinic acid
    rosma_e = Metabolite(
        'rosma_e',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='e')

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    amp_c = k12.metabolites.get_by_id("amp_c")

    #Diphosphate
    ppi_c = k12.metabolites.get_by_id("ppi_c")

    #CoA
    coa_c = k12.metabolites.get_by_id("coa_c")

    ############### Reactions ###############

    print("adding reactions")

    #caffeic acid transport reaction
    CAFFAt = Reaction('CAFFAt')
    CAFFAt.name = 'caffa transport'
    CAFFAt.subsystem = 'RA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        caffeic_acid_e: 1.0
    })

    CAFFCOA = uni.reactions.get_by_id("CAFFCOA")
    CAFFCOA.gene_reaction_rule = '(4CL)' #synthetic 4CL gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })


    #salvianic acid A + caffeoyl CoA -> Rosmarinic acid
    RAS = Reaction('RAS')
    RAS.name = 'rosmarinic acid synthase'
    RAS.subsystem = 'RA module'
    RAS.lower_bound = -1000  # This is the default
    RAS.upper_bound = 1000  # This is the default
    #Caffeoyl-CoA + 3-(3,4-Dihydroxyphenyl)lactate <=> CoA + Rosmarinate
    RAS.add_metabolites({
        caffcoa_c: -1.0,
        saa_c: -1.0,
        coa_c: 1.0,
        rosma_c: 1.0
    })
    RAS.gene_reaction_rule = '(ras)' #synthetic RAS gene


    #RA transport reaction
    RAt = Reaction('RAt')
    RAt.name = 'caffa transport'
    RAt.subsystem = 'RA module'
    RAt.lower_bound = -1000  # This is the default
    RAt.upper_bound = 1000  # This is the default
    RAt.add_metabolites({
        rosma_c: -1.0,
        rosma_e: 1.0
    })

    # add reactions to model
    MAM2.add_reactions([CAFFAt, CAFFCOA, SAAt, RAS, RAt])

    #exhange reactions for SA, caffa, and RA

    #MAM2.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'

    MAM2.add_boundary(MAM2.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    MAM2.add_boundary(MAM2.metabolites.get_by_id("saa_e"), type="exchange")
    MAM2.add_boundary(MAM2.metabolites.get_by_id("rosma_e"), type="exchange")

    #so that cell cannot eat rosmarinic acid
    MAM2.reactions.EX_rosma_e.lower_bound = 0

    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(MAM2), "Energy generation cycle detected"
    assert production_possible(MAM2, "EX_rosma_e"), "RA production not possible"

    #updtate the model ID
    MAM2.id = 'MAM2'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "MAM2.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(MAM2, file_path)


def sal9():

    print("loading GEMs...")

    print("... bl21 model")
    # read model for E. coli BL21
    bl21 = read_sbml_model('../GEMs/iHK1487.xml')

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    SAL9 = bl21.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #4-hydroxyphenyllactate
    hpl34_c = uni.metabolites.get_by_id("34hpl_c")
    hpl34_c.compartment = 'c'

    #3,4-dihydroxyphenylpyruvate
    dhpp34_c = Metabolite(
        '34dhpp_c',
        formula='C9H7O5',
        name='3,4-Dihydroxyphenylpyruvate',
        compartment='c')

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #loading existing metabolites for bl21

    #HPP
    hpp34_c = bl21.metabolites.get_by_id("34hpp_c")

    #NADH
    nadh_c = bl21.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = bl21.metabolites.get_by_id("nad_c")

    #H+
    h_c = bl21.metabolites.get_by_id("h_c")

    #H20
    h2o_c = bl21.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = bl21.metabolites.get_by_id("o2_c")

    ############### Reactions ###############

    print("adding reactions")

    #DLDH = uni.reactions.get_by_id("34HPLFM") #this reaction is for the mitochondiral compartment, written out for cytosol below

    #HPP -> 4-hydroxyphenyllactate
    DLDH = Reaction('DLDH')
    DLDH.name = 'D-lactate dehydrogenase'
    DLDH.subsystem = 'SAA module'
    DLDH.lower_bound = -1000  # This is the default
    DLDH.upper_bound = 1000  # This is the default
    #3,4-Hydroxyphenyllactate + NAD+ <=> 3,4-Hydroxyphenylpyruvate + NADH + H+
    DLDH.add_metabolites({
        hpl34_c: -1.0,
        nad_c: -1.0,
        hpp34_c: 1.0,
        nadh_c: 1.0,
        h_c: 1.0
    })
    DLDH.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #HPP -> 3,4-dihydroxyphenylpyruvate
    HPPHD = Reaction('HPPHD')
    HPPHD.name = '4-hydroxyphenylacetate 3-hydroxylase'
    HPPHD.subsystem = 'SAA module'
    HPPHD.lower_bound = -1000  # This is the default
    HPPHD.upper_bound = 1000  # This is the default
    #4-hydroxyphenylpyruvate + Oxygen + NADH + H+ <=> 3,4-dihydroxyphenylpyruvate + NAD+ + H2O
    HPPHD.add_metabolites({
        hpp34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        dhpp34_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPPHD.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaBC gene
    
    #3,4-dihydroxyphenylpyruvate -> salvianic acid A
    DHPPSA = Reaction('DHPPSA')
    DHPPSA.name = '4-hydroxyphenyllactate:NAD+ oxidoreductase'
    DHPPSA.subsystem = 'SAA module'
    DHPPSA.lower_bound = -1000  # This is the default
    DHPPSA.upper_bound = 1000  # This is the default
    #3,4-Dihydroxyphenyllactate + NAD+ <=> 3,4-Dihydroxyphenylpyruvate + NADH + H+
    DHPPSA.add_metabolites({
        saa_c: 1.0,
        nad_c: 1.0,
        dhpp34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0
    })
    DHPPSA.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #4-hydroxyphenyllactate -> salvianic acid A
    HPLSA = Reaction('HPLSA')
    HPLSA.name = '(R)-3-(4-hydroxyphenyl)lactate:NADP+ oxidoreductase'
    HPLSA.subsystem = 'SAA module'
    HPLSA.lower_bound = -1000  # This is the default
    HPLSA.upper_bound = 1000  # This is the default
    #4-Hydroxyphenyllactate + Oxygen + NADH + H+ <=> 3,4-Dihydroxyphenylacetic acid + NAD+ + H2O
    HPLSA.add_metabolites({
        hpl34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        saa_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPLSA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic ldh gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })

    # add reactions to model
    SAL9.add_reactions([DLDH, HPPHD, DHPPSA, HPLSA, SAAt])

    #exhange reactions for SA
    SAL9.add_boundary(SAL9.metabolites.get_by_id("saa_e"), type="exchange")

    #so that the cell cannot eat SA

    SAL9.reactions.EX_saa_e.lower_bound = 0
    
    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(SAL9), "Energy generation cycle detected"
    assert production_possible(SAL9, "EX_saa_e"), "SAA production not possible"

    #updtate the model ID
    SAL9.id = 'SAL9'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "SAL9.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(SAL9, file_path)


def cal2():

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model('../GEMs/iML1515.xml')

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    CAL2 = k12.copy()

    ############### Knock-outs ###############

    # pheA | gene b2599
    knock_out_model_genes(CAL2, ["b2599"])

    ############### Metabolites ###############

    print("adding metabolites")

    #p-Coumaric acid
    p_coumaric_acid_c = uni.metabolites.get_by_id("T4hcinnm_c")
    p_coumaric_acid_c.compartment = "c"

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = "c"

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = "e"

    #ammonia
    nh4_c = k12.metabolites.get_by_id("nh4_c")

    #tyrosine 
    tyr__L_c = k12.metabolites.get_by_id("tyr__L_c")

    #NADH
    nadh_c = k12.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = k12.metabolites.get_by_id("nad_c")

    #H+
    h_c = k12.metabolites.get_by_id("h_c")

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    adp_c = k12.metabolites.get_by_id("adp_c")

    #Diphosphate
    pi_c = k12.metabolites.get_by_id("pi_c")

    #H20
    h2o_c = k12.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = k12.metabolites.get_by_id("o2_c")

    ############### Reactions ###############

    print("adding reactions")

    #tyr -> p-coumaric acid
    TAL = Reaction('TAL')
    TAL.name = 'L-tyrosine ammonia-lyase'
    TAL.subsystem = 'CA module'
    TAL.lower_bound = -1000  # This is the default
    TAL.upper_bound = 1000  # This is the default
    #L-Tyrosine <=> 4-Coumarate + Ammonia
    TAL.add_metabolites({
        tyr__L_c: -1.0,
        p_coumaric_acid_c: 1.0,
        nh4_c: 1.0
    })
    TAL.gene_reaction_rule = 'tal' #synthetic rgTAL gene

    #p-coumaric acid -> caffeic acid
    COURCA = Reaction('COURCA')
    COURCA.name = '4-hydroxyphenylacetate 3-hydroxylase'
    COURCA.subsystem = 'CA module'
    COURCA.lower_bound = -1000  # This is the default
    COURCA.upper_bound = 1000  # This is the default
    #4-Coumarate <=> Caffeate
    COURCA.add_metabolites({
        p_coumaric_acid_c: -1.0,
        nadh_c: -1.0,
        o2_c: -1.0,
        h_c: -1.0,
        caffeic_acid_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    COURCA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaC gene

    #caffeic acid transport reaction
    CAFFAt = Reaction('34DHCINMt')
    CAFFAt.name = 'caffeic acid active transport'
    CAFFAt.subsystem = 'CA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        atp_c: -1.0,
        caffeic_acid_e: 1.0,
        adp_c: 1.0,
        pi_c: 1.0
    })

    # add reactions to model
    CAL2.add_reactions([TAL, COURCA, CAFFAt])

    #exhange reactions for CA
    #CAL2.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'
    CAL2.add_boundary(CAL2.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    
    #stopping caffeic acid uptake so that the model does not use this as substrate
    CAL2.reactions.EX_34dhcinm_e.lower_bound = 0
    
    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(CAL2), "Energy generation cycle detected"
    assert production_possible(CAL2, "EX_34dhcinm_e"), "CA production not possible"

    #updtate the model ID
    CAL2.id = 'CAL2'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "CAL2.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(CAL2, file_path)
    print("done")


def mam1():

    print("loading GEMs...")

    print("... bl21 model")
    # read model for E. coli BL21
    bl21 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "iHK1487.xml"))

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    MAM1 = bl21.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = 'c'

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = 'e'

    #caffeoyl CoA
    caffcoa_c = uni.metabolites.get_by_id("caffcoa_c")
    caffcoa_c.compartment = 'c'

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H10O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H10O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #rosmarinic acid
    rosma_c = Metabolite(
        'rosma_c',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='c')

    #external rosmarinic acid
    rosma_e = Metabolite(
        'rosma_e',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='e')

    #ATP
    atp_c = bl21.metabolites.get_by_id("atp_c")

    #AMP
    amp_c = bl21.metabolites.get_by_id("amp_c")

    #Diphosphate
    ppi_c = bl21.metabolites.get_by_id("ppi_c")

    #CoA
    coa_c = bl21.metabolites.get_by_id("coa_c")

    ############### Reactions ###############

    print("adding reactions")

    #caffeic acid transport reaction
    CAFFAt = Reaction('CAFFAt')
    CAFFAt.name = 'caffa transport'
    CAFFAt.subsystem = 'RA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        caffeic_acid_e: 1.0
    })

    CAFFCOA = uni.reactions.get_by_id("CAFFCOA")


    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })


    #salvianic acid A + caffeoyl CoA -> Rosmarinic acid
    RAS = Reaction('RAS')
    RAS.name = 'rosmarinic acid synthase'
    RAS.subsystem = 'RA module'
    RAS.lower_bound = -1000  # This is the default
    RAS.upper_bound = 1000  # This is the default
    #Caffeoyl-CoA + 3-(3,4-Dihydroxyphenyl)lactate <=> CoA + Rosmarinate
    RAS.add_metabolites({
        caffcoa_c: -1.0,
        saa_c: -1.0,
        coa_c: 1.0,
        rosma_c: 1.0
    })
    #RAS.gene_reaction_rule = '(ras)' #synthetic 4CL gene


    #RA transport reaction
    RAt = Reaction('RAt')
    RAt.name = 'caffa transport'
    RAt.subsystem = 'RA module'
    RAt.lower_bound = -1000  # This is the default
    RAt.upper_bound = 1000  # This is the default
    RAt.add_metabolites({
        rosma_c: -1.0,
        rosma_e: 1.0
    })

    # add reactions to model
    MAM1.add_reactions([CAFFAt, CAFFCOA, SAAt, RAS, RAt])

    #exhange reactions for SA, caffa, and RA

    #MAM1.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'

    MAM1.add_boundary(MAM1.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    MAM1.add_boundary(MAM1.metabolites.get_by_id("saa_e"), type="exchange")
    MAM1.add_boundary(MAM1.metabolites.get_by_id("rosma_e"), type="exchange")

    #so that cell cannot eat rosmarinic acid
    MAM1.reactions.EX_rosma_e.lower_bound = 0

    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(MAM1), "Energy generation cycle detected"
    assert production_possible(MAM1, "EX_rosma_e"), "RA production not possible"

    #updtate the model ID
    MAM1.id = 'MAM1'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "GEMs" / "MAM1.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(MAM1, file_path)


def mra():
    """Model for the MRA construct, that was grown in monoculture to express entire pathway"""

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model('../GEMs/iML1515.xml')

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    MRA = k12.copy()

    ############### Knock-outs ###############

    # pheA | gene b2599
    knock_out_model_genes(MRA, ["b2599"])

    ############### Metabolites ###############

    print("adding metabolites")

    #p-Coumaric acid
    p_coumaric_acid_c = uni.metabolites.get_by_id("T4hcinnm_c")
    p_coumaric_acid_c.compartment = "c"

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = "c"

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = "e"

    #ammonia
    nh4_c = k12.metabolites.get_by_id("nh4_c")

    #tyrosine 
    tyr__L_c = k12.metabolites.get_by_id("tyr__L_c")

    #NADH
    nadh_c = k12.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = k12.metabolites.get_by_id("nad_c")

    #H+
    h_c = k12.metabolites.get_by_id("h_c")

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    adp_c = k12.metabolites.get_by_id("adp_c")

    #Diphosphate
    pi_c = k12.metabolites.get_by_id("pi_c")

    #4-hydroxyphenyllactate
    hpl34_c = uni.metabolites.get_by_id("34hpl_c")
    hpl34_c.compartment = 'c'

    #3,4-dihydroxyphenylpyruvate
    dhpp34_c = Metabolite(
        '34dhpp_c',
        formula='C9H8O5',
        name='3,4-Dihydroxyphenylpyruvate',
        compartment='c')

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H10O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #HPP #TODO: check that this actually exists here
    hpp34_c = k12.metabolites.get_by_id("34hpp_c")

    #caffeoyl CoA
    caffcoa_c = uni.metabolites.get_by_id("caffcoa_c")
    caffcoa_c.compartment = 'c'

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H10O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #rosmarinic acid
    rosma_c = Metabolite(
        'rosma_c',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='c')

    #external rosmarinic acid
    rosma_e = Metabolite(
        'rosma_e',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='e')
    
    #CoA
    coa_c = k12.metabolites.get_by_id("coa_c")

    #H20
    h2o_c = k12.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = k12.metabolites.get_by_id("o2_c")

    ############### Reactions ###############

    print("adding reactions")

    #tyr -> p-coumaric acid
    TAL = Reaction('TAL')
    TAL.name = 'L-tyrosine ammonia-lyase'
    TAL.subsystem = 'CA module'
    TAL.lower_bound = -1000  # This is the default
    TAL.upper_bound = 1000  # This is the default
    #L-Tyrosine <=> 4-Coumarate + Ammonia
    TAL.add_metabolites({
        tyr__L_c: -1.0,
        p_coumaric_acid_c: 1.0,
        nh4_c: 1.0
    })
    #TAL.gene_reaction_rule = 'tal' #synthetic rgTAL gene

    #p-coumaric acid -> caffeic acid
    COURCA = Reaction('COURCA')
    COURCA.name = '4-hydroxyphenylacetate 3-hydroxylase'
    COURCA.subsystem = 'CA module'
    COURCA.lower_bound = -1000  # This is the default
    COURCA.upper_bound = 1000  # This is the default
    #4-Coumarate <=> Caffeate
    COURCA.add_metabolites({
        p_coumaric_acid_c: -1.0,
        nadh_c: -1.0,
        o2_c: -1.0,
        h_c: -1.0,
        caffeic_acid_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    COURCA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaC gene

    #HPP -> 4-hydroxyphenyllactate
    DLDH = Reaction('DLDH')
    DLDH.name = 'D-lactate dehydrogenase'
    DLDH.subsystem = 'SAA module'
    DLDH.lower_bound = -1000  # This is the default
    DLDH.upper_bound = 1000  # This is the default
    #3,4-Hydroxyphenyllactate + NAD+ <=> 3,4-Hydroxyphenylpyruvate + NADH + H+
    DLDH.add_metabolites({
        hpl34_c: -1.0,
        nad_c: -1.0,
        hpp34_c: 1.0,
        nadh_c: 1.0,
        h_c: 1.0
    })
    #DLDH.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #HPP -> 3,4-dihydroxyphenylpyruvate
    HPPHD = Reaction('HPPHD')
    HPPHD.name = '4-hydroxyphenylacetate 3-hydroxylase'
    HPPHD.subsystem = 'SAA module'
    HPPHD.lower_bound = -1000  # This is the default
    HPPHD.upper_bound = 1000  # This is the default
    #4-hydroxyphenylpyruvate + Oxygen + NADH + H+ <=> 3,4-dihydroxyphenylpyruvate + NAD+ + H2O
    HPPHD.add_metabolites({
        # hpp34_c: -1.0,
        # o2_c: -1.0,
        # nadh_c: -1.0,
        # h_c: -1.0,
        # dhpp34_c: 1.0,
        # nad_c: 1.0,
        # h2o_c: 1.0
        hpp34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        dhpp34_c: 1.0,
        nad_c: 1.0,
    })
    #HPPHD.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaBC gene
    
    #3,4-dihydroxyphenylpyruvate -> salvianic acid A
    DHPPSA = Reaction('DHPPSA')
    DHPPSA.name = '4-hydroxyphenyllactate:NAD+ oxidoreductase'
    DHPPSA.subsystem = 'SAA module'
    DHPPSA.lower_bound = -1000  # This is the default
    DHPPSA.upper_bound = 1000  # This is the default
    #3,4-Dihydroxyphenyllactate + NAD+ <=> 3,4-Dihydroxyphenylpyruvate + NADH + H+
    DHPPSA.add_metabolites({
        saa_c: 1.0,
        nad_c: 1.0,
        dhpp34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0
    })
    #DHPPSA.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #4-hydroxyphenyllactate -> salvianic acid A
    HPLSA = Reaction('HPLSA')
    HPLSA.name = '(R)-3-(4-hydroxyphenyl)lactate:NADP+ oxidoreductase'
    HPLSA.subsystem = 'SAA module'
    HPLSA.lower_bound = -1000  # This is the default
    HPLSA.upper_bound = 1000  # This is the default
    #4-Hydroxyphenyllactate + Oxygen + NADH + H+ <=> 3,4-Dihydroxyphenylacetic acid + NAD+ + H2O
    HPLSA.add_metabolites({
        # hpl34_c: -1.0,
        # o2_c: -1.0,
        # nadh_c: -1.0,
        # h_c: -1.0,
        # saa_c: 1.0,
        # nad_c: 1.0,
        # h2o_c: 1.0
        hpl34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        saa_c: 1.0,
        nad_c: 1.0,
    })
    #HPLSA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic ldh gene

    CAFFCOA = uni.reactions.get_by_id("CAFFCOA")

    #salvianic acid A + caffeoyl CoA -> Rosmarinic acid
    RAS = Reaction('RAS')
    RAS.name = 'rosmarinic acid synthase'
    RAS.subsystem = 'RA module'
    RAS.lower_bound = -1000  # This is the default
    RAS.upper_bound = 1000  # This is the default
    #Caffeoyl-CoA + 3-(3,4-Dihydroxyphenyl)lactate <=> CoA + Rosmarinate
    RAS.add_metabolites({
        caffcoa_c: -1.0,
        saa_c: -1.0,
        coa_c: 1.0,
        rosma_c: 1.0
    })
    #RAS.gene_reaction_rule = '(ras)' #synthetic 4CL gene


    #RA transport reaction
    RAt = Reaction('RAt')
    RAt.name = 'caffa transport'
    RAt.subsystem = 'RA module'
    RAt.lower_bound = -1000  # This is the default
    RAt.upper_bound = 1000  # This is the default
    RAt.add_metabolites({
        rosma_c: -1.0,
        rosma_e: 1.0
    })

    # add reactions to model
    MRA.add_reactions([TAL, COURCA, DLDH, HPPHD, DHPPSA, HPLSA, CAFFCOA, RAS, RAt])

    MRA.add_boundary(MRA.metabolites.get_by_id("rosma_e"), type="exchange")

    #so that cell cannot eat rosmarinic acid
    MRA.reactions.EX_rosma_e.lower_bound = 0

    print("Testing model")

    # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(MRA), "Energy generation cycle detected"
    assert production_possible(MRA, "EX_rosma_e"), "RA production not possible"

    #updtate the model ID
    MRA.id = 'MRA'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "GEMs" / "MRA.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(MRA, file_path)


def cal11():
    """GEM for CA module with knocked out glucose transport"""

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model('../GEMs/iML1515.xml')

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    CAL11 = k12.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #p-Coumaric acid
    p_coumaric_acid_c = uni.metabolites.get_by_id("T4hcinnm_c")
    p_coumaric_acid_c.compartment = "c"

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = "c"

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = "e"

    #ammonia
    nh4_c = k12.metabolites.get_by_id("nh4_c")

    #tyrosine 
    tyr__L_c = k12.metabolites.get_by_id("tyr__L_c")

    #NADH
    nadh_c = k12.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = k12.metabolites.get_by_id("nad_c")

    #H+
    h_c = k12.metabolites.get_by_id("h_c")

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    adp_c = k12.metabolites.get_by_id("adp_c")

    #Diphosphate
    pi_c = k12.metabolites.get_by_id("pi_c")

    #H20
    h2o_c = k12.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = k12.metabolites.get_by_id("o2_c")
    ############### Reactions ###############

    print("adding reactions")

    #tyr -> p-coumaric acid
    TAL = Reaction('TAL')
    TAL.name = 'L-tyrosine ammonia-lyase'
    TAL.subsystem = 'CA module'
    TAL.lower_bound = -1000  # This is the default
    TAL.upper_bound = 1000  # This is the default
    #L-Tyrosine <=> 4-Coumarate + Ammonia
    TAL.add_metabolites({
        tyr__L_c: -1.0,
        p_coumaric_acid_c: 1.0,
        nh4_c: 1.0
    })
    TAL.gene_reaction_rule = 'tal' #synthetic rgTAL gene

    #p-coumaric acid -> caffeic acid
    COURCA = Reaction('COURCA')
    COURCA.name = '4-hydroxyphenylacetate 3-hydroxylase'
    COURCA.subsystem = 'CA module'
    COURCA.lower_bound = -1000  # This is the default
    COURCA.upper_bound = 1000  # This is the default
    #4-Coumarate <=> Caffeate
    COURCA.add_metabolites({
        p_coumaric_acid_c: -1.0,
        nadh_c: -1.0,
        o2_c: -1.0,
        h_c: -1.0,
        caffeic_acid_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    COURCA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaC gene

    #caffeic acid transport reaction
    CAFFAt = Reaction('34DHCINMt')
    CAFFAt.name = 'caffeic acid active transport'
    CAFFAt.subsystem = 'CA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        atp_c: -1.0,
        caffeic_acid_e: 1.0,
        adp_c: 1.0,
        pi_c: 1.0
    })

    # add reactions to model
    CAL11.add_reactions([TAL, COURCA, CAFFAt])

    #exhange reactions for CA
    #CAL11.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'
    CAL11.add_boundary(CAL11.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    
    #stopping caffeic acid uptake so that the model does not use this as substrate
    CAL11.reactions.EX_34dhcinm_e.lower_bound = 0

    # allow uptake of xylose
    CAL11.reactions.EX_xyl__D_e.lower_bound = -10

    ############### Knock-outs ###############

    # glucose transport reactions
    CAL11.reactions.GLCtex_copy1.bounds = (0.0, 0.0)
    CAL11.reactions.GLCtex_copy2.bounds = (0.0, 0.0)

    # # pheA | gene b2599
    # knock_out_model_genes(CAL11, ["b2599"])
    # # ptsH | b2415
    # knock_out_model_genes(CAL11, ["b2415"])
    # # ptsI | b2416
    # knock_out_model_genes(CAL11, ["b2416"])
    # # crr | b2417
    # knock_out_model_genes(CAL11, ["b2417"])
    # # ydiB | b1692
    # knock_out_model_genes(CAL11, ["b1692"])
    
    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(CAL11), "Energy generation cycle detected"
    assert production_possible(CAL11, "EX_34dhcinm_e"), "CA production not possible"

    #updtate the model ID
    CAL11.id = 'CAL11'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "CAL11.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(CAL11, file_path)


def mam3():
    """GEM for RA module with knocked out glucose transport"""

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "iML1515.xml"))

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    MAM3 = k12.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = 'c'

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = 'e'

    #caffeoyl CoA
    caffcoa_c = uni.metabolites.get_by_id("caffcoa_c")
    caffcoa_c.compartment = 'c'

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #rosmarinic acid
    rosma_c = Metabolite(
        'rosma_c',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='c')

    #external rosmarinic acid
    rosma_e = Metabolite(
        'rosma_e',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='e')

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    amp_c = k12.metabolites.get_by_id("amp_c")

    #Diphosphate
    ppi_c = k12.metabolites.get_by_id("ppi_c")

    #CoA
    coa_c = k12.metabolites.get_by_id("coa_c")

    ############### Reactions ###############

    print("adding reactions")

    #caffeic acid transport reaction
    CAFFAt = Reaction('CAFFAt')
    CAFFAt.name = 'caffa transport'
    CAFFAt.subsystem = 'RA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        caffeic_acid_e: 1.0
    })

    # 4CL reaction CA -> caffeoyl-CoA
    CAFFCOA = uni.reactions.get_by_id("CAFFCOA")
    CAFFCOA.gene_reaction_rule = '(4CL)' #synthetic 4CL gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })


    #salvianic acid A + caffeoyl CoA -> Rosmarinic acid
    RAS = Reaction('RAS')
    RAS.name = 'rosmarinic acid synthase'
    RAS.subsystem = 'RA module'
    RAS.lower_bound = -1000  # This is the default
    RAS.upper_bound = 1000  # This is the default
    #Caffeoyl-CoA + 3-(3,4-Dihydroxyphenyl)lactate <=> CoA + Rosmarinate
    RAS.add_metabolites({
        caffcoa_c: -1.0,
        saa_c: -1.0,
        coa_c: 1.0,
        rosma_c: 1.0
    })
    RAS.gene_reaction_rule = '(ras)' #synthetic moRAS gene


    #RA transport reaction
    RAt = Reaction('RAt')
    RAt.name = 'caffa transport'
    RAt.subsystem = 'RA module'
    RAt.lower_bound = -1000  # This is the default
    RAt.upper_bound = 1000  # This is the default
    RAt.add_metabolites({
        rosma_c: -1.0,
        rosma_e: 1.0
    })

    # add reactions to model
    MAM3.add_reactions([CAFFAt, CAFFCOA, SAAt, RAS, RAt])

    #exhange reactions for SA, caffa, and RA

    #MAM3.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'

    MAM3.add_boundary(MAM3.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    MAM3.add_boundary(MAM3.metabolites.get_by_id("saa_e"), type="exchange")
    MAM3.add_boundary(MAM3.metabolites.get_by_id("rosma_e"), type="exchange")

    #so that cell cannot eat rosmarinic acid
    MAM3.reactions.EX_rosma_e.lower_bound = 0

    # allow uptake of xylose
    MAM3.reactions.EX_xyl__D_e.lower_bound = -10

    ############### Knock-outs ###############

    # glucose transport reactions
    MAM3.reactions.GLCtex_copy1.bounds = (0.0, 0.0)
    MAM3.reactions.GLCtex_copy2.bounds = (0.0, 0.0)

    # # pheA | gene b2599
    # knock_out_model_genes(CAL11, ["b2599"])
    # # ptsH | b2415
    # knock_out_model_genes(CAL11, ["b2415"])
    # # ptsI | b2416
    # knock_out_model_genes(CAL11, ["b2416"])
    # # crr | b2417
    # knock_out_model_genes(CAL11, ["b2417"])
    # # ydiB | b1692
    # knock_out_model_genes(CAL11, ["b1692"])

    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(MAM3), "Energy generation cycle detected"
    assert production_possible(MAM3, "EX_rosma_e"), "RA production not possible"

    #updtate the model ID
    MAM3.id = 'MAM3'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "MAM3.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(MAM3, file_path)


def sal11():
    """ GEM for SAA module with knocked out xylose transport"""

    print("loading GEMs...")

    print("... bl21 model")
    # read model for E. coli BL21
    bl21 = read_sbml_model('../GEMs/iHK1487.xml')

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    SAL11 = bl21.copy()

    ############### Metabolites ###############

    print("adding metabolites")

    #4-hydroxyphenyllactate
    hpl34_c = uni.metabolites.get_by_id("34hpl_c")
    hpl34_c.compartment = 'c'

    #3,4-dihydroxyphenylpyruvate
    dhpp34_c = Metabolite(
        '34dhpp_c',
        formula='C9H7O5',
        name='3,4-Dihydroxyphenylpyruvate',
        compartment='c')

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #loading existing metabolites for bl21

    #HPP
    hpp34_c = bl21.metabolites.get_by_id("34hpp_c")

    #NADH
    nadh_c = bl21.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = bl21.metabolites.get_by_id("nad_c")

    #H+
    h_c = bl21.metabolites.get_by_id("h_c")

    #H20
    h2o_c = bl21.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = bl21.metabolites.get_by_id("o2_c")

    ############### Reactions ###############

    print("adding reactions")

    #DLDH = uni.reactions.get_by_id("34HPLFM") #this reaction is for the mitochondiral compartment, written out for cytosol below

    #HPP -> 4-hydroxyphenyllactate
    DLDH = Reaction('DLDH')
    DLDH.name = 'D-lactate dehydrogenase'
    DLDH.subsystem = 'SAA module'
    DLDH.lower_bound = -1000  # This is the default
    DLDH.upper_bound = 1000  # This is the default
    #3,4-Hydroxyphenyllactate + NAD+ <=> 3,4-Hydroxyphenylpyruvate + NADH + H+
    DLDH.add_metabolites({
        hpl34_c: -1.0,
        nad_c: -1.0,
        hpp34_c: 1.0,
        nadh_c: 1.0,
        h_c: 1.0
    })
    DLDH.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #HPP -> 3,4-dihydroxyphenylpyruvate
    HPPHD = Reaction('HPPHD')
    HPPHD.name = '4-hydroxyphenylacetate 3-hydroxylase'
    HPPHD.subsystem = 'SAA module'
    HPPHD.lower_bound = -1000  # This is the default
    HPPHD.upper_bound = 1000  # This is the default
    #4-hydroxyphenylpyruvate + Oxygen + NADH + H+ <=> 3,4-dihydroxyphenylpyruvate + NAD+ + H2O
    HPPHD.add_metabolites({
        hpp34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        dhpp34_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPPHD.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaBC gene
    
    #3,4-dihydroxyphenylpyruvate -> salvianic acid A
    DHPPSA = Reaction('DHPPSA')
    DHPPSA.name = '4-hydroxyphenyllactate:NAD+ oxidoreductase'
    DHPPSA.subsystem = 'SAA module'
    DHPPSA.lower_bound = -1000  # This is the default
    DHPPSA.upper_bound = 1000  # This is the default
    #3,4-Dihydroxyphenyllactate + NAD+ <=> 3,4-Dihydroxyphenylpyruvate + NADH + H+
    DHPPSA.add_metabolites({
        saa_c: 1.0,
        nad_c: 1.0,
        dhpp34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0
    })
    DHPPSA.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #4-hydroxyphenyllactate -> salvianic acid A
    HPLSA = Reaction('HPLSA')
    HPLSA.name = '(R)-3-(4-hydroxyphenyl)lactate:NADP+ oxidoreductase'
    HPLSA.subsystem = 'SAA module'
    HPLSA.lower_bound = -1000  # This is the default
    HPLSA.upper_bound = 1000  # This is the default
    #4-Hydroxyphenyllactate + Oxygen + NADH + H+ <=> 3,4-Dihydroxyphenylacetic acid + NAD+ + H2O
    HPLSA.add_metabolites({
        hpl34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        saa_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPLSA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaBC gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })

    # add reactions to model
    SAL11.add_reactions([DLDH, HPPHD, DHPPSA, HPLSA, SAAt])

    # exhange reactions for SA
    SAL11.add_boundary(SAL11.metabolites.get_by_id("saa_e"), type="exchange")

    # ############### Knock-outs ###############

    # xylose uptake
    SAL11.reactions.XYLtex.bounds = (0.0, 0.0)

    # # xylA (b3565) | ECD_03417
    # knock_out_model_genes(SAL11, ["ECD_03417"])

    # so that the cell cannot eat SA
    SAL11.reactions.EX_saa_e.lower_bound = 0

    
    print("Testing model")

    # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(SAL11), "Energy generation cycle detected"
    assert production_possible(SAL11, "EX_saa_e"), "SAA production not possible"

    # updtate the model ID
    SAL11.id = 'SAL11'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "SAL11.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(SAL11, file_path)
    

def rau2():
    """E. coli K12 PH2 carrying the CA module for co-culture."""

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "iML1515.xml"))

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    RAU2 = k12.copy()

    ############### Knock-outs ###############

    # pheA | gene b2599
    knock_out_model_genes(RAU2, ["b2599"])

    ############### Metabolites ###############

    print("adding metabolites")

    #p-Coumaric acid
    p_coumaric_acid_c = uni.metabolites.get_by_id("T4hcinnm_c")
    p_coumaric_acid_c.compartment = "c"

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = "c"

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = "e"

    #ammonia
    nh4_c = k12.metabolites.get_by_id("nh4_c")

    #tyrosine 
    tyr__L_c = k12.metabolites.get_by_id("tyr__L_c")

    #NADH
    nadh_c = k12.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = k12.metabolites.get_by_id("nad_c")

    #H+
    h_c = k12.metabolites.get_by_id("h_c")

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    adp_c = k12.metabolites.get_by_id("adp_c")

    #Diphosphate
    pi_c = k12.metabolites.get_by_id("pi_c")

    #H20
    h2o_c = k12.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = k12.metabolites.get_by_id("o2_c")

    ############### Reactions ###############

    print("adding reactions")

    #tyr -> p-coumaric acid
    TAL = Reaction('TAL')
    TAL.name = 'L-tyrosine ammonia-lyase'
    TAL.subsystem = 'CA module'
    TAL.lower_bound = -1000  # This is the default
    TAL.upper_bound = 1000  # This is the default
    #L-Tyrosine <=> 4-Coumarate + Ammonia
    TAL.add_metabolites({
        tyr__L_c: -1.0,
        p_coumaric_acid_c: 1.0,
        nh4_c: 1.0
    })
    TAL.gene_reaction_rule = 'tal' #synthetic rgTAL gene

    #p-coumaric acid -> caffeic acid
    COURCA = Reaction('COURCA')
    COURCA.name = '4-hydroxyphenylacetate 3-hydroxylase'
    COURCA.subsystem = 'CA module'
    COURCA.lower_bound = -1000  # This is the default
    COURCA.upper_bound = 1000  # This is the default
    #4-Coumarate <=> Caffeate
    COURCA.add_metabolites({
        p_coumaric_acid_c: -1.0,
        nadh_c: -1.0,
        o2_c: -1.0,
        h_c: -1.0,
        caffeic_acid_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    COURCA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaC gene

    #caffeic acid transport reaction
    CAFFAt = Reaction('34DHCINMt')
    CAFFAt.name = 'caffeic acid active transport'
    CAFFAt.subsystem = 'CA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        atp_c: -1.0,
        caffeic_acid_e: 1.0,
        adp_c: 1.0,
        pi_c: 1.0
    })

    # add reactions to model
    RAU2.add_reactions([TAL, COURCA, CAFFAt])

    #exhange reactions for CA
    #RAU2.metabolites.get_by_id("34dhcinm_e").compartment = 'e' #set compartment to be 'e' rather than 'C_e'
    RAU2.add_boundary(RAU2.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    
    #stopping caffeic acid uptake so that the model does not use this as substrate
    RAU2.reactions.EX_34dhcinm_e.lower_bound = 0
    
    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(RAU2), "Energy generation cycle detected"
    assert production_possible(RAU2, "EX_34dhcinm_e"), "CA production not possible"

    #updtate the model ID
    RAU2.id = 'RAU2'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "RAU2.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(RAU2, file_path)


def rad4():
    """E. coli K12 PH2 carrying the SAA and RA module for co-culture."""

    print("loading GEMs...")

    print("... k12 model")
    # read model for E. coli K12
    k12 = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "iML1515.xml"))

    print("... universal model")
    #read the universal model
    uni = read_sbml_model(str(ROOT_DIR / "community_modelling" / "GEMs" / "bigg_universe.xml"))

    RAD4 = k12.copy()

    ############### Knock-outs ###############

    # pheA | gene b2599
    knock_out_model_genes(RAD4, ["b2599"])

    # tyrB | gene b4054
    knock_out_model_genes(RAD4, ["b4054"])

    ############### Metabolites ###############

    print("adding metabolites")

    #4-hydroxyphenyllactate
    hpl34_c = uni.metabolites.get_by_id("34hpl_c")
    hpl34_c.compartment = 'c'

    #caffeic acid
    caffeic_acid_c = uni.metabolites.get_by_id("34dhcinm_c")
    caffeic_acid_c.compartment = 'c'

    #external caffeic acid
    caffeic_acid_e = uni.metabolites.get_by_id("34dhcinm_e")
    caffeic_acid_e.compartment = 'e'

    #caffeoyl CoA
    caffcoa_c = uni.metabolites.get_by_id("caffcoa_c")
    caffcoa_c.compartment = 'c'

    #3,4-dihydroxyphenylpyruvate
    dhpp34_c = Metabolite(
        '34dhpp_c',
        formula='C9H7O5',
        name='3,4-Dihydroxyphenylpyruvate',
        compartment='c')

    #3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_c = Metabolite(
        'saa_c',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='c')

    #external 3,4-Dihydroxyphenyllactic acid | Salvianic acid A | SAA
    saa_e = Metabolite(
        'saa_e',
        formula='C9H9O5',
        name='3,4-Dihydroxyphenyllactic acid',
        compartment='e')

    #loading existing metabolites for K12 TODO: check that these actually exist in the k12 model!!

    #HPP
    hpp34_c = k12.metabolites.get_by_id("34hpp_c")

    #NADH
    nadh_c = k12.metabolites.get_by_id("nadh_c")

    #NAD+
    nad_c = k12.metabolites.get_by_id("nad_c")

    #H+
    h_c = k12.metabolites.get_by_id("h_c")

    #H20
    h2o_c = k12.metabolites.get_by_id("h2o_c")

    #O2
    o2_c = k12.metabolites.get_by_id("o2_c")

    #rosmarinic acid
    rosma_c = Metabolite(
        'rosma_c',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='c')

    #external rosmarinic acid
    rosma_e = Metabolite(
        'rosma_e',
        formula='C18H16O8',
        name='Rosmarinic acid',
        compartment='e')

    #ATP
    atp_c = k12.metabolites.get_by_id("atp_c")

    #AMP
    amp_c = k12.metabolites.get_by_id("amp_c")

    #Diphosphate
    ppi_c = k12.metabolites.get_by_id("ppi_c")

    #CoA
    coa_c = k12.metabolites.get_by_id("coa_c")

    ############### Reactions ###############

    print("adding reactions")

    #HPP -> 4-hydroxyphenyllactate
    DLDH = Reaction('DLDH')
    DLDH.name = 'D-lactate dehydrogenase'
    DLDH.subsystem = 'SAA module'
    DLDH.lower_bound = -1000  # This is the default
    DLDH.upper_bound = 1000  # This is the default
    #3,4-Hydroxyphenyllactate + NAD+ <=> 3,4-Hydroxyphenylpyruvate + NADH + H+
    DLDH.add_metabolites({
        hpl34_c: -1.0,
        nad_c: -1.0,
        hpp34_c: 1.0,
        nadh_c: 1.0,
        h_c: 1.0
    })
    DLDH.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #HPP -> 3,4-dihydroxyphenylpyruvate
    HPPHD = Reaction('HPPHD')
    HPPHD.name = '4-hydroxyphenylacetate 3-hydroxylase'
    HPPHD.subsystem = 'SAA module'
    HPPHD.lower_bound = -1000  # This is the default
    HPPHD.upper_bound = 1000  # This is the default
    #4-hydroxyphenylpyruvate + Oxygen + NADH + H+ <=> 3,4-dihydroxyphenylpyruvate + NAD+ + H2O
    HPPHD.add_metabolites({
        hpp34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        dhpp34_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPPHD.gene_reaction_rule = '(hpaB and hpaC)' #synthetic hpaBC gene
    
    #3,4-dihydroxyphenylpyruvate -> salvianic acid A
    DHPPSA = Reaction('DHPPSA')
    DHPPSA.name = '4-hydroxyphenyllactate:NAD+ oxidoreductase'
    DHPPSA.subsystem = 'SAA module'
    DHPPSA.lower_bound = -1000  # This is the default
    DHPPSA.upper_bound = 1000  # This is the default
    #3,4-Dihydroxyphenyllactate + NAD+ <=> 3,4-Dihydroxyphenylpyruvate + NADH + H+
    DHPPSA.add_metabolites({
        saa_c: 1.0,
        nad_c: 1.0,
        dhpp34_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0
    })
    DHPPSA.gene_reaction_rule = 'ldh' #synthetic ldh gene

    #4-hydroxyphenyllactate -> salvianic acid A
    HPLSA = Reaction('HPLSA')
    HPLSA.name = '(R)-3-(4-hydroxyphenyl)lactate:NADP+ oxidoreductase'
    HPLSA.subsystem = 'SAA module'
    HPLSA.lower_bound = -1000  # This is the default
    HPLSA.upper_bound = 1000  # This is the default
    #4-Hydroxyphenyllactate + Oxygen + NADH + H+ <=> 3,4-Dihydroxyphenylacetic acid + NAD+ + H2O
    HPLSA.add_metabolites({
        hpl34_c: -1.0,
        o2_c: -1.0,
        nadh_c: -1.0,
        h_c: -1.0,
        saa_c: 1.0,
        nad_c: 1.0,
        h2o_c: 1.0
    })
    HPLSA.gene_reaction_rule = '(hpaB and hpaC)' #synthetic ldh gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })

    #caffeic acid transport reaction
    CAFFAt = Reaction('CAFFAt')
    CAFFAt.name = 'caffa transport'
    CAFFAt.subsystem = 'RA module'
    CAFFAt.lower_bound = -1000  # This is the default
    CAFFAt.upper_bound = 1000  # This is the default
    CAFFAt.add_metabolites({
        caffeic_acid_c: -1.0,
        caffeic_acid_e: 1.0
    })

    # 4CL reaction CA -> caffeoyl-CoA
    CAFFCOA = uni.reactions.get_by_id("CAFFCOA")
    CAFFCOA.gene_reaction_rule = '(4CL)' #synthetic 4CL gene

    #SAA transport reaction
    SAAt = Reaction('SAAt')
    SAAt.name = 'caffa transport'
    SAAt.subsystem = 'RA module'
    SAAt.lower_bound = -1000  # This is the default
    SAAt.upper_bound = 1000  # This is the default
    SAAt.add_metabolites({
        saa_c: -1.0,
        saa_e: 1.0
    })

    #salvianic acid A + caffeoyl CoA -> Rosmarinic acid
    RAS = Reaction('RAS')
    RAS.name = 'rosmarinic acid synthase'
    RAS.subsystem = 'RA module'
    RAS.lower_bound = -1000  # This is the default
    RAS.upper_bound = 1000  # This is the default
    #Caffeoyl-CoA + 3-(3,4-Dihydroxyphenyl)lactate <=> CoA + Rosmarinate
    RAS.add_metabolites({
        caffcoa_c: -1.0,
        saa_c: -1.0,
        coa_c: 1.0,
        rosma_c: 1.0
    })
    RAS.gene_reaction_rule = '(ras)' #synthetic moRAS gene


    #RA transport reaction
    RAt = Reaction('RAt')
    RAt.name = 'caffa transport'
    RAt.subsystem = 'RA module'
    RAt.lower_bound = -1000  # This is the default
    RAt.upper_bound = 1000  # This is the default
    RAt.add_metabolites({
        rosma_c: -1.0,
        rosma_e: 1.0
    })

    RAD4.add_reactions([DLDH, HPPHD, DHPPSA, HPLSA, SAAt, CAFFAt, CAFFCOA, RAS, RAt])

    # exhange reactions for SA, CA, RA
    RAD4.add_boundary(RAD4.metabolites.get_by_id("34dhcinm_e"), type="exchange")
    RAD4.add_boundary(RAD4.metabolites.get_by_id("saa_e"), type="exchange")
    RAD4.add_boundary(RAD4.metabolites.get_by_id("rosma_e"), type="exchange")

    #so that cell cannot eat rosmarinic acid
    RAD4.reactions.EX_rosma_e.lower_bound = 0

    print("Testing model")

    # # run tests to see that the model functions as it should
    assert no_energy_gen_cycles(RAD4), "Energy generation cycle detected"
    assert production_possible(RAD4, "EX_rosma_e"), "RA production not possible"

    #updtate the model ID
    RAD4.id = 'RAD4'

    # write model to file
    file_path = ROOT_DIR / "community_modelling" / "RAsynthesis" / "GEMs" / "RAD4.xml"
    print("Writing model to file as", str(file_path))
    write_sbml_model(RAD4, file_path)


if __name__ == '__main__':
    freeze_support()
    #mam2()
    #sal9()
    #cal2()
    #cal11()
    #mam3()
    #sal11()
    #rau2()
    rad4()
