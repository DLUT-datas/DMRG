
# construct weighted differential network
import networkx as nx
import pickle
from itertools import islice 


# read gene symbol from dataset
filename = 'gene_id_protein.txt' # 
gene_symbol=[]    # gene symbol
protein_id_info=[]      # protein ID

with open(filename)as f:
    for line in islice(f,1,None):
        line=line.strip().split('\t')
        gene_symbol.append(line[0])
        protein_id_info.append(line[1])
# pickle.dump(gene_symbol, 
#             open( 'gene_symbol.pickle', "wb"))

# read PPI network
g_ppi_all = nx.read_graphml("g_ppi_all.graphml")

# combine PPI network node info and differential network node info
# Unified node info
n = len(gene_symbol)
sim_score_ppi =[]
edge_info_ppi = []
for i in range(n):
    temp_id_1 = protein_id_info[i]
    print(i)
    if temp_id_1 != '0':
        for j in range(n):
            if i>=j:
                continue
            else:
                temp_id_2 = protein_id_info[i]
                if temp_id_2 != '0':
                    
                    if g_ppi_all.has_edge(protein_id_info[i],protein_id_info[j]):

                        a = g_ppi_all.get_edge_data(protein_id_info[i],protein_id_info[j])
                        sim_score_ppi.append(a['weight'])
                        edge_info_ppi.append([gene_symbol[i],gene_symbol[j]])

import pickle

pickle.dump(sim_score_ppi, 
            open( 'sim_score_ppi.pickle', "wb"))
pickle.dump(edge_info_ppi, 
            open( 'edge_info_ppi.pickle', "wb"))


    
# PPI differential network graph
import networkx as nx
g_ppi_diff=nx.Graph()
# add nodes and edges
g_ppi_diff.add_nodes_from(gene_symbol)
g_ppi_diff.add_edges_from(edge_info_ppi)
# add edge weights from PPI score
num = len(sim_score_ppi)
for i in range(num):
    g_ppi_diff.add_edge(edge_info_ppi[i][0], edge_info_ppi[i][1], weight=sim_score_ppi[i])




# read differential correlation score

import pickle
outfile =  'sim_score_cor_diff.pickle'
with open(outfile, 'rb') as file:
    diff_cor_score = pickle.load(file)

# compute the edge weigths
num = len(diff_cor_score)
summary_score = []
for i in range(num):
    if abs(diff_cor_score[i])<0.8:
        summary_score.append(0)
    else:
        summary_score.append(abs(diff_cor_score[i])*sim_score_ppi[i])

# weighted differential network

g_weighted_diff=nx.Graph()
# add nodes and edges 
g_weighted_diff.add_nodes_from(gene_symbol)
g_weighted_diff.add_edges_from(edge_info_ppi)
# add edge weight

num = len(sim_score_ppi)
for i in range(num):
    g_weighted_diff.add_edge(edge_info_ppi[i][0], edge_info_ppi[i][1], weight=summary_score[i])

nx.write_graphml_lxml(g_weighted_diff, "g_weighted_diff_network.graphml")
