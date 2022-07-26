# -*- coding: utf-8 -*-

# output information

import networkx as nx
import csv
import pickle
import scipy.io as scio

# read node info
outfile =  'gene_symbol.pickle'
with open(outfile, 'rb') as file:
    gene_symbol = pickle.load(file)

# read network graph
g_weighted_diff = nx.read_graphml("g_weighted_diff_network.graphml")

g_sum_node_degree = g_weighted_diff.degree(weight='weight')

# output top 10 node info and neighbors

# select top-rank 10 nodes based on node weighted degree
node_degree_value = []
for one_node in gene_symbol:
    node_degree_value.append(g_sum_node_degree[one_node])
    
sort_degree = sorted(enumerate(node_degree_value),key=lambda x:x[1], reverse=True)
sort_degree_id = [m[0] for m in sort_degree]
sort_degree_value = [m[1] for m in sort_degree]

net_node = [gene_symbol[i] for i in sort_degree_id[0:10] ]

# net_node = ['MYC','CD44','CDK6','WNT5A','ALB','EPHA4','BLNK','PTEN','RANBP2','MTOR']

#net_node = ['TRIM25']
#RANBP2
for i in net_node:
    out_name = i +'g_sum_network_hub_node_info.txt'
    
    output_info = open(out_name,'a',newline='',encoding='utf-8')    
    csv_write = csv.writer(output_info,dialect='excel')
   
    edlist=[]
    for ed in g_weighted_diff.edges():
         if i in ed:
             temp_link = g_weighted_diff.get_edge_data(ed[0],ed[1])
             if temp_link['weight'] !=0:
                 if ed[0] not in edlist:
                     edlist.append(ed[0])
                     csv_write.writerow([ed[0]])
                 if ed[1] not in edlist:
                     edlist.append(ed[1])
                     csv_write.writerow([ed[1]])
                   
              
    output_info.close() 