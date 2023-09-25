"""
Replicating the plots from the paper.
For dynamic sim results.
"""

import seaborn as sns
import pandas as pd

def plot_relative_abundance(biomass_df: pd.DataFrame):

    total_BM = biomass_df["CAL11"] + biomass_df["SAL11"] + biomass_df["MAM3"]
    CAL11_frac = biomass_df["CAL11"] / total_BM
    SAL11_frac = biomass_df["SAL11"] / total_BM
    MAM3_frac = biomass_df["MAM3"] / total_BM

    df = pd.concat([total_BM, CAL11_frac, SAL11_frac, MAM3_frac], axis=1, keys=["total_BM", "CAL11", "SAL11", "MAM3"])
    df["cycles"] = df.index
    plot_df = df.melt(id_vars="cycles", value_vars=["CAL11", "SAL11", "MAM3"], value_name="relative_abundance", var_name="strain")

    sns.lineplot(plot_df, x="cycles", y="relative_abundance", hue="strain", hue_order=["CAL11", "SAL11", "MAM3"])
    plt.ylabel("Subpopulation fraction")
    return plt