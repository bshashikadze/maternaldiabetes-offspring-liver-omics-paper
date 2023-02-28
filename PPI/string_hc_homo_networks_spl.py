import numpy as np
import networkx as nx
import itertools as it
import random as rd
import pickle as pk
import pandas as pd
from collections import (defaultdict,Counter)

def network_spl(network):
    network_lcc = network.subgraph(max(nx.connected_components(network), key=len))
    spl={}
    pairwise=it.combinations(network_lcc.nodes(), 2)
    for pair in pairwise:
        spl[pair]= nx.shortest_path_length(network_lcc, pair[0], pair[1])
    return spl

string_id_symbol_converter={}

with open("9606.protein.info.v11.5.txt") as f:
    next(f)
    for line in f:
        string_id_symbol_converter[line.split("\t")[0][5:]]=line.split("\t")[1]

homo_string_full_dict={}
with open("9606.protein.links.v11.5.txt") as f:
    next(f)
    for line in f:
        homo_string_full_dict[line.split(" ")[0][5:],line.split(" ")[1][5:]]=float(line.split(" ")[2].strip())


homo_string_full_dict_symbol={}
for pp, c_score in homo_string_full_dict.items():
    homo_string_full_dict_symbol[string_id_symbol_converter[pp[0]],string_id_symbol_converter[pp[1]]]=c_score


G_homo_string_full=nx.Graph()
G_homo_string_full_hc=nx.Graph()

for pp, c_score in homo_string_full_dict_symbol.items():
    G_homo_string_full.add_edge(pp[0],pp[1])

for pp, c_score in homo_string_full_dict_symbol.items():
    if c_score>700:
        G_homo_string_full_hc.add_edge(pp[0],pp[1])

print(G_homo_string_full.number_of_nodes())
print(G_homo_string_full_hc.number_of_nodes())


networkspl=network_spl(G_homo_string_full_hc)
with open('output/spl_G_homo_string_full_hc_fast.pickle', 'wb') as handle:
    pk.dump(networkspl, handle, protocol=pk.HIGHEST_PROTOCOL)
