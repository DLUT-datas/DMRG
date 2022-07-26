# -*- coding: utf-8 -*-

# construct PPI network 


from itertools import islice 
import numpy as np
import csv
import pickle
import networkx as nx

def find_outlier_value(protein_score):

#    max_value = max(protein_score)
#    ave_value = np.mean(protein_score)
#    std_value = np.std(protein_score)
    down_precentie_value = np.percentile(protein_score,25)
#    median_value = np.median(protein_score)
    up_precentie_value = np.percentile(protein_score,75)

    Q_range = up_precentie_value-down_precentie_value
    down_outlier = down_precentie_value + 1.5*Q_range
    up_outlier = up_precentie_value + 1.5*Q_range
    
    return down_outlier,up_outlier

def write_txt(file_name,unique_node):
    f = open(file_name,'a', newline='',encoding = 'utf-8')

    for data in unique_node:
    
        f.write(data+'\n')
        
    f.close()

def discretization(protein_score):
    #<Q1,Q1-Q3,Q3-up_outlier,up_outlier>
#    num = len(protein_score)
    disc_protein_score = []
    Q1 = np.percentile(protein_score,25)
    Q3 = np.percentile(protein_score,75)
    Q_range = Q3-Q1
    up_outlier = Q3 + 1.5*Q_range
    for score in protein_score:
        if score<Q1:
            disc_protein_score.append(0.25)
        elif score>=Q1 and  score<Q3:
            disc_protein_score.append(0.5)
        elif score>=Q3 and  score<up_outlier:
            disc_protein_score.append(0.75)
        else:
            disc_protein_score.append(1)
    return disc_protein_score

# read protein info
filename = '9606.protein.info.v11.0 (1).txt'
protein_id_name = []
protein_id = []
protein_name = []
with open(filename) as f:
    for line in islice(f,1,None):
        line=line.strip().split('\t')
        protein_id.append(line[0])
        protein_name.append(line[1])
        protein_id_name.append([line[0],line[1]])
        
# read PPI network info 
filename = '9606.protein.links.v11.0 (1).txt' # 
protein_1=[]    # protein ID
protein_2=[]      # gene symbol
protein_score=[]      # score
protein_link=[]
with open(filename)as f:
    for line in islice(f,1,None):
        line=line.strip().split(' ')
        protein_1.append(line[0])
        protein_2.append(line[1])
        protein_link.append([line[0],line[1]])
        protein_score.append(float(line[2]))
       
# score normalization[ 0.25,0.5,0.75,1]
norm_protein_score = discretization(protein_score)

# network node
node = protein_1+protein_2
unique_node = list(set(node))
unique_node_info = []
for i in range(len(unique_node)):
    loc = protein_id.index(unique_node[i])
#    print(protein_id.index(unique_node[0]))
    unique_node_info.append([protein_id[loc],protein_name[loc]])

#pickle.dump(unique_node_info, 
#            open( 'ppi_unique_node_info.pickle', "wb"))
#pickle.dump(unique_node, 
#           open( 'ppi_unique_node.pickle', "wb"))   
write_txt('unique_node.txt',unique_node)
len(unique_node)

# PPI network graph
g_ppi_all=nx.Graph()# null graph
# add node and edge
g_ppi_all.add_nodes_from(unique_node)
g_ppi_all.add_edges_from(protein_link)

# add weight (score)
num = len(norm_protein_score)
for i in range(num):
    g_ppi_all.add_edge(protein_1[i], protein_2[i], weight=norm_protein_score[i])

nx.write_graphml_lxml(g_ppi_all, "g_ppi_all.graphml")




