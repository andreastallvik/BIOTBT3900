"""Kinetic parameter dictionaries."""

KINETIC_PARAMS = {
    "M5": {
        "km": {
            "EX_xylan4_e": 0.01,
            "EX_xylan8_e": 0.01,
            "EX_xyl__D_e": 1,
            'EX_val__L_e': 0.5,
            'EX_arg__L_e': 0.5,
            'EX_asp__L_e': 0.5,
            'EX_dhptd_e': 0.5,
            'EX_glu__L_e': 0.5,
            'EX_ile__L_e': 0.5,
            'EX_ser__L_e': 0.5,
            'EX_thr__L_e': 0.5,
            'EX_ala__L_e': 0.5,
            'EX_cys__L_e': 0.5,
            'EX_gly_e': 0.5,
            'EX_his__L_e': 0.5,
            'EX_leu__L_e': 0.5,
            'EX_met__L_e': 0.5,
            'EX_phe__L_e': 0.5,
            'EX_pro__L_e': 0.5,
            'EX_tyr__L_e': 0.5,
            'EX_trp__L_e': 0.5,
            'EX_lys__L_e': 0.5,
        },
        "vmax": {
            "EX_xylan4_e": 10,
            "EX_xylan8_e": 6,
            "EX_xyl__D_e": 6,
        },
        "hill": {
            "EX_xylan4_e": 1,
            "EX_xylan8_e": 1,
            "EX_xyl__D_e": 1,
        },
    },
    "NJ4": {
        "km": {
            "EX_ac_e": 0.01,
            "EX_but_e": 0.01,
            "EX_xyl__D_e": 1,
            'EX_val__L_e': 0.5,
            'EX_arg__L_e': 0.5,
            'EX_asp__L_e': 0.5,
            'EX_dhptd_e': 0.5,
            'EX_glu__L_e': 0.5,
            'EX_ile__L_e': 0.5,
            'EX_ser__L_e': 0.5,
            'EX_thr__L_e': 0.5,
            'EX_ala__L_e': 0.5,
            'EX_cys__L_e': 0.5,
            'EX_gly_e': 0.5,
            'EX_his__L_e': 0.5,
            'EX_leu__L_e': 0.5,
            'EX_met__L_e': 0.5,
            'EX_phe__L_e': 0.5,
            'EX_pro__L_e': 0.5,
            'EX_tyr__L_e': 0.5,
            'EX_trp__L_e': 0.5,
            'EX_lys__L_e': 0.5,
        },
        "vmax": {
            "EX_ac_e": 10,
            "EX_but_e": 10,
            "EX_xyl__D_e": 6,
        },
        "hill": {
            "EX_xylan4_e": 1,
            "EX_xylan8_e": 1,
            "EX_xyl__D_e": 1,
        },
    },
}

NJ4_VMAX = KINETIC_PARAMS["NJ4"]["vmax"]
NJ4_KM = KINETIC_PARAMS["NJ4"]["km"]
NJ4_HILL = KINETIC_PARAMS["NJ4"]["hill"]
M5_VMAX = KINETIC_PARAMS["M5"]["vmax"]
M5_KM = KINETIC_PARAMS["M5"]["km"]
M5_HILL = KINETIC_PARAMS["M5"]["hill"]
