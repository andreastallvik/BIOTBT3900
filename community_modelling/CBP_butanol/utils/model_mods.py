def make_ABE_pathway_irreverable(model):
    """Makes the reactions of the ABE pathway irreversible to adhere to the observed phenotype."""

    irr_model = model.copy()
    
    reactions = ["POR_syn",
                "ACACT1r",
                "HACD1",
                "ECOAH1",
                "ACOAD1fr",
                "ACOAD1",
                "BTCOARx",
                "PBUTT",
                "ADCi",
                "PTAr"]
    
    reverse_reactions = ["ALCD4", "BUTKr", "BUTCT2", "ACKr", "ACACCT", "ACALD"]

    for rx in reactions:
        model.reactions.get_by_id(rx).bounds = (0, 1000)
    
    for rx in reverse_reactions:
        model.reactions.get_by_id(rx).bounds = (-1000, 0)

    return irr_model