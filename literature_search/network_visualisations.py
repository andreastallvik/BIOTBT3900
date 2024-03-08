"""
TODO: create a more compact bipartite network using the substrate map and product map
"""

import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import itertools as it
from collections import Counter
from pyvis.network import Network


def organism_co_occurance(data: pd.DataFrame, filename: str, genus_level: bool = True):
    """Creates co-occurance graph for organisms, either at species of genus (default) level. Saves network to html."""

    df = data.copy()

    #read data from file and crate datastructures in lists
    nodes = []
    edges = []

    #iterate over the rows of the dataframe
    for index, row in df.iterrows():

        #get lists of organisms
        organisms = row["Organisms"].split(", ")

        if genus_level:
            #taking the first word to get the genus instead of the organisms
            genuses = [organism.split(' ')[0] for organism in organisms]
            organisms = genuses

        #get a list of each edge
        co_occurrences = list(it.combinations(organisms, 2))

        #add all the organisms the the node list
        nodes.extend(organisms)

        #add all the egdes to the global edge list
        edges.extend(co_occurrences)

        #add a self-loop if there is only a single organism listed
        if len(organisms) == 1:
            edges.append((organisms[0], organisms[0]))

    #build network
    coOccurrenceNetwork = nx.Graph()

    #add nodes w/ size
    for node in nodes:
        s = 5 + nodes.count(node)
        coOccurrenceNetwork.add_node(node, size = s)

    #add edges w/weights
    edge_counts = Counter(edges)
    for e in edge_counts:
        u = e[0]
        v = e[1]
        c = edge_counts[e]
        coOccurrenceNetwork.add_edge(u, v, weight = c)

    #visualise network
    nt = Network('1000px', '1800px') #, select_menu=True)
    nt.from_nx(coOccurrenceNetwork)
    filepath = "plots/" + filename + ".html"
    nt.show(filepath, notebook=False)

def bipartite_organism_compounds(data: pd.DataFrame, filename: str, genus_level: bool = True):

    df = data.copy()

    #read data from file and crate datastructures in lists
    compound_nodes = []
    organism_nodes = []
    edges = []


    #iterate over the rows of the dataframe
    for index, row in df.iterrows():

        #get lists of substrates, products, and organisms
        substrates = row["Substrate"].split(", ")
        products = row["Product"].split(", ")
        organisms = row["Organisms"].split(", ")

        if genus_level:
            #taking the first word to get the genus instead of the organisms
            genuses = [organism.split(' ')[0] for organism in organisms]
            organisms = genuses

        for substrate in substrates:
            compound_nodes.append(substrate)

            for organism in organisms:
                organism_nodes.append(organism)
                #add edge from the substrate to the organism
                edges.append((substrate, organism))
                
                for product in products:
                    compound_nodes.append(product)
                    #add edge from the organism to the product
                    edges.append((organism, product))
                 

    #build network
    bipartiteNetwork = nx.DiGraph()

    for node in compound_nodes:
        s = 5 + compound_nodes.count(node)
        bipartiteNetwork.add_node(node, bipartite = 0, color = "red", size = s)

    for node in organism_nodes:
        s = 5 + organism_nodes.count(node)
        bipartiteNetwork.add_node(node, bipartite = 1, shape = 'square', size = s)

    bipartiteNetwork.add_edges_from(edges)

    #visualise network
    nt = Network('1000px', '1800px', directed=True)
    nt.from_nx(bipartiteNetwork)
    filepath = "plots/" + filename + ".html"
    nt.show(filepath, notebook=False)