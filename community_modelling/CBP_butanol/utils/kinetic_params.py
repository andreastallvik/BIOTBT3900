"""Kinetic parameter dictionaries."""

#TODO: replace each instance of "None" with correct values
# CoAT params from that 1999 article (but v. approximate)

KINETIC_PARAMS = {
    "M5": {
        "km": {
            "EX_xylan4_e": 5,
            "EX_xylan8_e": 5,
            "EX_xyl__D_e": 5,
        },
        "vmax": {
            "EX_xylan4_e": 10,
            "EX_xylan8_e": 10,
            "EX_xyl__D_e": 10,
        }
    },
    "NJ4": {
        "km": {
            "EX_ac_e": 60,
            "EX_but_e": 161,
            "EX_xyl__D_e": 5,
        },
        "vmax": {
            "EX_ac_e": 10,
            "EX_but_e": 10,
            "EX_xyl__D_e": 10,
        }
    },
}

NJ4_VMAX = KINETIC_PARAMS["NJ4"]["vmax"]
NJ4_KM = KINETIC_PARAMS["NJ4"]["km"]
M5_VMAX = KINETIC_PARAMS["M5"]["vmax"]
M5_KM = KINETIC_PARAMS["M5"]["km"]
