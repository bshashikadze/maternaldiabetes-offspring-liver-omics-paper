import numpy as np
import networkx as nx
import itertools as it
import random as rd
import pickle as pk
import pandas as pd


def random_selection(G,size):
    while True:
        rnd_sample=rd.sample(list(G.nodes()),size)
        if 'nan' in rnd_sample:
            rnd_sample=rd.sample(list(G.nodes()),size)
        else:
            break
    return rnd_sample

def calculate_closest_distance(G_ppi_lcc,spl, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            if node_from==node_to:
                val =0
            else:
                try:
                    try:
                        val = spl[node_from,node_to]
                    except:
                        val = spl[node_to,node_from]
                except:
                    val=len(nx.shortest_path(G_ppi_lcc,source=node_from, target=node_to))
            values.append(val)
        d = min(values)
        #print d,
        values_outer.append(d)
    d = np.mean(values_outer)
    #print d
    return d

def calculate_proximity(network,spl, nodes_from, nodes_to, exp_random):   #let's calculate the proximity between the two sets
    d=calculate_closest_distance(network,spl,nodes_from,nodes_to)         #exposure and disease genesets
    m=exp_random[len(nodes_from),len(nodes_to)][0]
    s=exp_random[len(nodes_from),len(nodes_to)][1]
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s) #(z, pval)

#Here, we define the pig networks

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
G_homo_string_full_hc_lcc = G_homo_string_full_hc.subgraph(max(nx.connected_components(G_homo_string_full_hc), key=len))

#Here, we define the set of degs for each condition

with open('output/homo_proximal_geneset_dict.pickle', 'rb') as handle:
    homo_proximal_geneset_dict = pk.load(handle)

#Here we define the set of diseases

with open('input/disgenet_gda_curated_score.pickle', 'rb') as handle:
    disgenet_gda_curated_score = pk.load(handle)

disgenet_disease_associations_df=pd.read_csv("input/disease_associations.tsv",sep="\t")
disgenet_gda_curated_score_name={}

for dis_id,geneset in disgenet_gda_curated_score.items():
    disgenet_gda_curated_score_name[disgenet_disease_associations_df[disgenet_disease_associations_df['diseaseId']==dis_id]['diseaseName'].values[0]]=geneset


#Let's create a new dictionary that consider for each netwotk the diseases that are present in that network

disgenet_gda_curated_score_multigenes_dict={}
for dis,geneset in disgenet_gda_curated_score_name.items():
    new_geneset=set(geneset) & G_homo_string_full_hc_lcc.nodes()
    if len(new_geneset)>1:    #at least two genes associated with the disease
        disgenet_gda_curated_score_multigenes_dict[dis]=new_geneset


#Let's import the pre-calculated spl dictionary
with open('output/spl_G_homo_string_full_hc_fast.pickle', 'rb') as handle:
    spl_G_homo_string_full_hc_fast = pk.load(handle)


#Let's create a disease geneset distribution for each network

disease_geneset_distribution=[]
for geneset in disgenet_gda_curated_score_multigenes_dict.values():
        disease_geneset_distribution.append(len(geneset))


homo_dis_proximity={}
for condition,geneset in homo_proximal_geneset_dict.items():
    homo_degs=set(geneset)&G_homo_string_full_hc_lcc.nodes()
    homo_spec_size=len(homo_degs)
    homo_random_distance={}

    for dis_size in set(disease_geneset_distribution):
        S = 10000
        l_rnd_dist = []
        for s in range(S):
            rnd_exp_set=random_selection(G_homo_string_full_hc_lcc,homo_spec_size)
            rnd_dis_set=random_selection(G_homo_string_full_hc_lcc,dis_size)
            rnd_dist=calculate_closest_distance(G_homo_string_full_hc_lcc,spl_G_homo_string_full_hc_fast,rnd_exp_set,rnd_dis_set)
            l_rnd_dist.append(rnd_dist)
        mu = np.mean(l_rnd_dist)
        std = np.std(l_rnd_dist)

        homo_random_distance[homo_spec_size,dis_size]=[mu,std]

    for dis in disgenet_gda_curated_score_multigenes_dict.keys():
        homo_dis_proximity[condition,dis]=calculate_proximity(G_homo_string_full_hc_lcc,spl_G_homo_string_full_hc_fast,homo_degs,disgenet_gda_curated_score_multigenes_dict[dis],homo_random_distance)

#Let's save it
with open('output/homo_degs_network_proximity.pickle', 'wb') as handle:
    pk.dump(homo_dis_proximity, handle, protocol=pk.HIGHEST_PROTOCOL)
