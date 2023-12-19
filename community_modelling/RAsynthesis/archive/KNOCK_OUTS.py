"""
Lists of knock-outs that need be performed.

for example in cobrapy:

from cobra.manipulation import knock_out_model_genes
import KNOCK_OUTS

knock_out_model_genes(SAL11, KNOCK_OUTS.SAL11_KO)
knock_out_model_genes(MAM3, KNOCK_OUTS.MAM3_KO)
knock_out_model_genes(CAL11, KNOCK_OUTS.CAL11_KO)

or in reframed:

import KNOCK_OUTS


"""

CAL2_KO = ["b2599"]

MRA_KO = ["b2599"]

CAL11_KO = ["b2599", "b2415", "b2416", "b2417", "b1692"]

MAM3_KO = ["b2599", "b2415", "b2416", "b2417", "b1692", "b3281"]

SAL11_KO = ["ECD_03417"]