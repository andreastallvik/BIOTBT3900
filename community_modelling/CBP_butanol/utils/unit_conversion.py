"""Unit conversion functions."""

# dictionary of molar mass for each metabolite
MOLAR_MASS = {"xyl__D_e": 150.13, 
        "etoh_e": 46.068, 
        "but_e": 88.11, 
        "btoh_e": 74.12, 
        "ac_e": 59.044, 
        "acetone_e": 58.08,
        "xylan4_e": 600.52,
        "xylan8_e": 1201.04,
        "butanol": 74.12,
        "acetate": 59.044,
        "ethanol": 46.068,
        "butyrate": 88.11,
        "acetone": 58.08,
        "xylose": 150.13}

def mmol_to_g_per_L(met_name, met_mmol, volume = 0.05):
    """Convert mmol to g/L. Can convert a single value or an array.""" 
    
    # divide by 1000 (-> mol), multiply by molar mass (-> g) and divide by volume (-> g/L)
    return (met_mmol / 1000) * MOLAR_MASS[met_name] / volume

def g_per_L_to_mmol(met_name, met_conc, volume = 0.05):
    """Convert g/L to mmol. Can convert a single value or an array.""" 

    # multiply by volume (-> g), divide by molar mass (-> mol) and multiply by 1000 (-> mmol)
    return (met_conc * volume) / MOLAR_MASS[met_name] * 1000
