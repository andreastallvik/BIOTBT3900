"""
Functions to generate the plots from literature search.

TODO: move all the plotting code (that I wanna keep) from defined_plots and sankey_diagram here, use notebook to survey results
in this notebook make short comments on the results from the figures, recap:)
maybe also network things here? -> nah use different 
"""

import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import random
from maps import SUBSTRATE_MAP, PRODUCT_MAP

#plt.style.use('seaborn-v0_8-pastel')

def publications_per_year(data: pd.DataFrame, save_fig = False):
    """Bar graph showing publications per year"""

    #TODO: add functionality to plot only a given range of years, look at how this was done in stacked_publications_hist

    df = data.copy()

    #add in the years without any publications
    publications_per_year = df.groupby("Year").size()
    all_years = pd.Series(0, index=range(1981,2023))
    publications_per_year = publications_per_year.add(all_years, fill_value=0)

    #make plot
    publications_per_year.plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/publications_per_year.svg", format="svg")


def stacked_publications_hist(data: pd.DataFrame, start_year = 2000, save_fig=False):
    """Stacked bar graph showing publications per year for defined / undefined communities"""

    df = data.copy()

    grouped = df.groupby(["Year", "Community type"])
    counts = grouped["Title"].count()

    idx = pd.MultiIndex.from_product([range(start_year,2023), ['Defined', 'Undefined']])

    # Use the reindex method to create a new series with the new index
    counts = counts.reindex(idx)

    # Use the fillna method to fill in any missing values with 0
    counts = counts.fillna(0)

    #make plot
    counts.unstack().plot(kind="bar", stacked=True)

    plt.style.use('default')
    plt.rcParams["figure.figsize"] = [6.5, 5]
    plt.rcParams["figure.autolayout"] = True

    plt.ylabel("Number of publications")
    plt.legend(title='Community type')
    plt.title("Publications using microbial comunities as production platforms")
    plt.show()

    if save_fig:
        plt.savefig("plots/stacked_publications_hist.svg", format="svg")


def publication_by_journal(data: pd.DataFrame, save_fig = False):
    """Bar chart showing publications per year"""

    df = data.copy()

    #Plot showing publications by journal
    publications_per_journal = df.groupby("Journal").size()

    #sort by decending order
    publications_per_journal = publications_per_journal.sort_values(ascending=False)

    #make the plot
    publications_per_journal.plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/publication_by_journal.svg", format="svg")


def substrate_hist(df: pd.DataFrame, save_fig = False):
    """Pareto histogram of substrates"""

    s_df = df.copy()

    s_df["Substrate"] = s_df["Substrate"].str.split(", ")

    # Use the explode function to create a new row for each substrate
    s_df = s_df.explode("Substrate")

    #use the substrate map to get more generalised results
    s_df["Substrate"] = [SUBSTRATE_MAP[s] for s in s_df["Substrate"]]

    #if having used the substrate map -> remove the duplicates for same papers
    for index, row in s_df.iterrows():
        all_substrates = row["Substrate"]
        unique_substrates = set(all_substrates)
        row["Substrate"] = list(unique_substrates)

    #group
    substrate_per_journal = s_df.groupby("Substrate").size()

    #sort by decending order
    substrate_per_journal = substrate_per_journal.sort_values(ascending=False)

    #make the plot
    substrate_per_journal.plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/substrate_hist.svg", format="svg")


def product_hist(df: pd.DataFrame, save_fig = False):
    """Pareto histogram of products"""

    p_df = df.copy()

    p_df["Product"] = p_df["Product"].str.split(", ")

    # Use the explode function to create a new row for each animal
    p_df = p_df.explode("Product")

    #use the product map to get more generalised results
    p_df["Product"] = [PRODUCT_MAP[p] for p in p_df["Product"]]

    #if having used the product map -> remove the duplicates for same papers
    for index, row in p_df.iterrows():
        all_products = row["Product"]
        unique_products = set(all_products)
        row["Product"] = list(unique_products)

    #group
    product_per_journal = p_df.groupby("Product").size()

    #sort by decending order
    product_per_journal = product_per_journal.sort_values(ascending=False)

    #make the plot
    product_per_journal.plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/product_hist.svg", format="svg")


def genus_hist(df: pd.DataFrame, save_fig = False):
    """Pareto histogram of organisms at genus level"""

    g_df = df.copy()

    g_df["Organisms"] = g_df["Organisms"].str.split(", ")

    #remove the duplicates for papers that use co-cultrues with the same genuses
    for index, row in g_df.iterrows():
        all_genuses = row["Organisms"]
        unique_genuses = set(all_genuses)
        row["Organisms"] = list(unique_genuses)

    # Use the explode function to create a new row for each animal
    g_df = g_df.explode("Organisms")

    #get only the genus
    g_df['Organisms'] = g_df['Organisms'].str.split().str[0]

    #group
    genus_per_journal = g_df.groupby("Organisms").size()

    #sort by decending order
    genus_per_journal = genus_per_journal.sort_values(ascending=False)

    #make the plot
    genus_per_journal[:13].plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/genus_hist.svg", format="svg")


def species_hist(df: pd.DataFrame, n: int, save_fig = False):
    """Pareto histogram of the n most common organisms """

    o_df = df.copy()

    o_df["Organisms"] = o_df["Organisms"].str.split(", ")

    # Use the explode function to create a new row for each animal
    o_df = o_df.explode("Organisms")

    #group
    organism_per_journal = o_df.groupby("Organisms").size()

    #sort by decending order
    organism_per_journal = organism_per_journal.sort_values(ascending=False)

    #make the plot
    organism_per_journal[:n].plot.bar()
    plt.ylabel("Number of publications")
    plt.show()

    if save_fig:
        plt.savefig("plots/species_hist.svg", format="svg")


def pick_colour(name: str, accent: str = None):
    """Helper function to pick colour for sankey diagram."""
    
    color_wheel = 'grey'

    if accent:
        color_wheel = px.colors.qualitative.Alphabet[8]

        if accent == "Clostridium":
            if name == "Clostridium":
                color_wheel = px.colors.qualitative.Pastel[6]

        if accent == "Escherichia":    
            if name == "Escherichia":
                color_wheel = px.colors.qualitative.Pastel[1]

    else:
        if name == "Clostridium":
                color_wheel = px.colors.qualitative.Pastel[6]
        elif name == "Escherichia":
                color_wheel = px.colors.qualitative.Pastel[1]
        else:
            color_wheel = random.choice(px.colors.qualitative.Pastel)

    return color_wheel


def sankey_genus(df: pd.DataFrame, color_accent = None, save_fig: bool = False, filename: str = None):
    """Sankey diagram showing substrate -> genus -> product.

    Args:
        df (pd.DataFrame): data
        color_accent (str, optional): "Escherichia" or "Clostridium", colours only the accent and leaves the rest gray. Defaults to None.
        save_fig (bool, optional): Save fig to pdf. Defaults to False.
        filename (str, optional): Filename if saving file. Defaults to None.

    """

    #NOTE: see sankey_diagram for more variations on sankey diagram plots.

    #creating a new dataframe to use in building the sankey chart

    #create an empty dataframe
    sankey_data = pd.DataFrame({'Source': [],
                    'Target': [],
                    'Value': [],
                    'Color': []})

    organism_colouring = dict()
    colour_counter = 0

    #iterate over the rows of the dataframe
    for index, row in df.iterrows():

        #get lists of substrates, products, and organisms
        all_substrates = row["Substrate"].split(", ")
        organisms = row["Organisms"].split(", ")
        all_products = row["Product"].split(", ")

        substrates = [SUBSTRATE_MAP[substrate] for substrate in all_substrates]
        products = [PRODUCT_MAP[product] for product in all_products]

        #taking the first word to get the genus instead of the organisms
        genuses = [organism.split(' ')[0] for organism in organisms]

        #adds colour to each organism using the pick_colour function
        for genus in genuses:
            if genus not in organism_colouring:
                organism_colouring[genus] = pick_colour(genus, color_accent)

        #add links from each substrate to each organism
        for substrate in substrates:
            for genus in genuses:
                sankey_data.loc[len(sankey_data.index)] = [substrate, genus, 1, organism_colouring[genus]]

        #add links from each organism to each substrate
        for genus in genuses:
            for product in products:
                sankey_data.loc[len(sankey_data.index)] = [genus, product, 1, organism_colouring[genus]]

    #needed modification for building the plot

    #get each unique source_target and a mapping to their index
    #unique_source_target = list(pd.unique(sankey_data[['Source', 'Target']].values.ravel('K')))
    unique_source_target = list(pd.unique(sankey_data[['Source', 'Target', 'Color']].values.ravel('K')))
    mapping_dict = {k: v for v, k in enumerate(unique_source_target)}
    sankey_data['Source'] = sankey_data['Source'].map(mapping_dict)
    sankey_data['Target'] = sankey_data['Target'].map(mapping_dict)
    sankey_dict = sankey_data.to_dict(orient='list')

    #setting colours for the nodes
    node_colours = [organism_colouring.get(i, "grey") for i in unique_source_target]

    #create diagram
    import plotly.express as px

    fig = go.Figure(data=[go.Sankey(
        orientation = "h",
        node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = unique_source_target,
        #color = "grey"
        color = node_colours
        ),
        link = dict(
        source = sankey_dict["Source"],
        target = sankey_dict["Target"],
        value = sankey_dict["Value"],
        color = sankey_dict["Color"]
        #color = [px.colors.qualitative.Plotly[unique_source_target.index(i) % len(px.colors.qualitative.Plotly)] for i in unique_source_target]
    ))])
    
    layout = dict(
        #title = "Sankey Diagram for consortia uses",
    height = 850,
    width = 1300,
    font = dict(
      size = 8),)

    #show figure
    fig.update_layout(layout)
    fig.show()

    if save_fig:
        if not filename:
            raise TypeError("Missing filename, cannot save figure.")
        
        filepath = "plots/" + filename + ".pdf"
        
        fig.write_image(filepath)
        fig.write_image(filepath)

        # NOTE: Saving the fig twice is on purpose! 
        # Due to a bug in plotly that puts a "Loading [MathJax]/extensions/MathMenu.js" box on the first plot, I am saving the fig twice, overwriting the first one
