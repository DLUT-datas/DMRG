
# construct differential network
import scipy.io as scio
import scipy.stats as stats
import numpy
from itertools import islice 
import csv
import numpy as np
import pickle
import networkx as nx

# read gene symbol and protein ID
filename = 'gene_id_protein.txt' # 
gene_symbol=[]    # gene symbol 
protein_id_info=[]      #已知蛋白质对应的基因名字

with open(filename)as f:
    for line in islice(f,1,None):
        line=line.strip().split('\t')
        gene_symbol.append(line[0])
        protein_id_info.append(line[1])


dataFile = 'prostate_tran.mat'
data = scio.loadmat(dataFile)


train_data = data['train_data']
train_label = data['train_label']
train_data_pos = train_data[0:24,:]
train_data_neg = train_data[24:,:]




m,n = train_data.shape
sim_score_pos = []
sim_score_neg = []
sim_score_link = []
sim_score_cor_diff = []
for i in range(n):
    temp_pos_1 = train_data_pos[:,i]
    temp_neg_1 = train_data_neg[:,i]
    print(i)
    for j in range(n):
        if i>=j:
            continue
        else:
            
            temp_pos_2 = train_data_pos[:,j]
            temp_neg_2 = train_data_neg[:,j]
            temp_sim_score_pos,p = stats.spearmanr(temp_pos_1, temp_pos_2)
            temp_sim_score_neg,p= stats.spearmanr(temp_neg_1, temp_neg_2)

            # sim_score_pos.append(temp_sim_score_pos)
            # sim_score_neg.append(temp_sim_score_neg)
            temp_diff_cor_score = abs(temp_sim_score_pos  - temp_sim_score_neg)
            sim_score_cor_diff.append(temp_diff_cor_score)
            if temp_diff_cor_score >0.8:
                sim_score_link.append([gene_symbol[i],gene_symbol[j],temp_diff_cor_score])

                
pickle.dump(sim_score_link, 
                open( 'sim_score_link.pickle', "wb"))

pickle.dump(sim_score_cor_diff, 
                open( 'sim_score_cor_diff.pickle', "wb"))

# differential network graph
g_diff_network =nx.Graph()# null graph
# add node and edge
g_diff_network.add_nodes_from(gene_symbol)
g_diff_network.add_edges_from(sim_score_link)

# add weight (score)
num = len(sim_score_link)
for ed in sim_score_link:
    g_diff_network.add_edge(ed[0], ed[1], weight=ed[2])

nx.write_graphml_lxml(g_diff_network, "g_diff_network.graphml")





