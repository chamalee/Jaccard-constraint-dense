
"""
Facebook dataset
"""

import networkx as nx
import time
from datetime import datetime, timedelta
from collections import OrderedDict
import matplotlib.pyplot as plt
import sys
import operator
import os
import copy
import numpy as np
import utils
import copy
import time
import os
import pandas as pd
from dsalgo import *
from find_densest_distinct_sets import *

filepath = os.path.join("Data","facebook.txt")
edgesTS, nodes, edges =  utils.readFile(filepath)

df = pd.DataFrame(edgesTS)
df.columns = ['source','target','timestamp']
# header_list =  ['source','target', 'timestamp']

# Remove null value
df = df[df['target'].isnull() != True]
df = df[df['source'].isnull() != True]
df = df[df['timestamp'].isnull() != True]
#sort
df=df.sort_values('timestamp')


def _swap (row):
    if row['source'] > row['target']:
        row['source'] , row['target'] =row['target'] , row['source']
    return row

# Undirected graph
df=df.apply(lambda row: _swap(row), axis=1)
#scale timestamps for zeroth reference point
refValue = df['timestamp'].min()
df['timestamp'] -= refValue

# Remove self-loops
df = df[((df['source'] ) != (df['target']))]
keys = df.timestamp.unique()
data_dic = {k: [] for k in keys}

print('TAU')
print(len(data_dic.keys()))

for _row in df.values:
        data_dic[_row[2]].append((_row[0],_row[1]))

snapshots =  []

sub = []
deg = []
nodes = set()

del_keys = []
def checkLimitedEdgeTimestamp(_limit):
    for i, (k, v) in enumerate(data_dic.items()):
        if len(set(v)) < _limit:
            del_keys.append(k)
              
_limit = 5
checkLimitedEdgeTimestamp(_limit)


for k in del_keys:
    del data_dic[k]

print('TAU')
print(len(data_dic.keys()))
#%% snaps
edg_lst = []
avg_nodes = 0
avg_edges = 0
for k, ls in data_dic.items():
    G=nx.Graph()  
    G.add_edges_from(ls, weight = 1)
    snapshots.append(G)
    
    # print(len(G.nodes()))
    avg_nodes += len(G.nodes())
    # print(len(G.edges()))
    avg_edges += len(G.edges()) 
        
    nodes = nodes.union( set(G.nodes()))
    edg_lst += [e for e in G.edges]
    
numberOfGraphs = len(snapshots)

avg_nodes /= numberOfGraphs 
avg_edges /= numberOfGraphs 

print('avg edges: %s'%avg_edges)
print('nodes:  {}'.format(len(nodes)))

#%%  algo
comb = snapshots[0]
for i in range(1, len(snapshots)):

    H = snapshots[i]
    comb = nx.compose(comb, H)

nodes = comb.nodes()

[obj1, subg1] =ip_based_dcs_sum(0, snapshots, comb)

densities = []
for i in range(len(snapshots)):
    G = snapshots[i]
    sub_g = G.subgraph(subg1)
    den = sub_g.number_of_edges()/len(subg1)
    densities.append(den)
    # print(den, end=", ")
print('IP-based')
print(sum(densities))
dcs = set(subg1)
dcs_den = sum(densities)
  #%%
class GraphObj:

    def __init__(self) :
        # Store the edge list as a dictionary (edge:weight)
        self.w_d = {}

    # Assuming that the edge is bidirectional
    def AddEdge (self, src , dst , weight ) :
    
        if (src,dst) in self.w_d:
              self.w_d[(src,dst)] += weight
        elif (dst, src) in self.w_d:
              self.w_d[(dst, src)] += weight
        else:
              self.w_d[(src,dst)] = weight        

    def get_wd (self):
        return  self.w_d

  #%% create avg graph

g = GraphObj()
for  x, y in edg_lst:
    g.AddEdge(x,y,1)
    
w_d =  g.get_wd()  
G_avg=nx.Graph()

for idx, val in w_d.items():
    G_avg.add_edge(idx[0],idx[1],weight=val)

 #%% Other methods
print('dcs density LP') 
d, induced , dcs1 = densest_subgraph_w(G_avg)
print(len(dcs1))
den = 0        
for key in range(0,numberOfGraphs):
    subgrpah_snap = snapshots[key].subgraph(dcs1)
    den+= (subgrpah_snap.number_of_edges()/ len(dcs1) )
print(den)

dcs2, d2 = exact_densest_weighted(G_avg, numberOfGraphs)
print('dcs density weighted') 
print(len(dcs2))
den = 0        
for key in range(0,numberOfGraphs):
    subgrpah_snap = snapshots[key].subgraph(dcs2)
    den+= (subgrpah_snap.number_of_edges()/ len(dcs2) )
print(den)

 #%%

ld = []
local_sum_den = 0
for snap in snapshots:
    d, induced ,  sol = densest_subgraph_w(snap)
    local_sum_den += d
    ld.append(set(sol))
print('local sum denisty')    
print(local_sum_den )
 #%%

arr = [0.3,0.5,0.7]

algo_ver = 0    
k =  len(snapshots)


for item  in arr:
    print('----------------------')
    
    if algo_ver == 0:
        # hard constraint
        alp = item
        set_dic =  hard(snapshots, dcs,  alp)  
    
    elif algo_ver == 1: 
        # soft contraint - 1

        print(item)
        lam  = item*dcs_den/k
        # lam  = 2*lam*den/(k*(k-1))        
        print(lam)
        
        [s1,set_dic1] = soft_1_1(snapshots, ld, lam)
        [s2, set_dic2] = soft_1(snapshots, dcs, lam)
        print('lam  {} s1: {} s2: {}'.format(lam,s1,s2))
        
        if s1 > s2 :
            set_dic = set_dic1
        else:
            set_dic = set_dic2
    elif algo_ver == 2: 
        # soft contraint - 2
        # lam = 5
        print(item)
        lam  = item*dcs_den/k
        print('lam ::: ',lam)
        set_dic = soft_2(snapshots, dcs, lam, nodes)
    
    else:
        print('no algo found')
    
    
 #    Results

    print('sets density') 
    den= 0     
    d_l = []      
    for key, val in set_dic.items():   
        # print('size {} set {}'.format(len(val),val))
        # print('size {} '.format(len(val)))
        den+=comDensity(snapshots[key],val)
        d_l.append(comDensity(snapshots[key],val))
    # print(d_l)
    # print(min(d_l))
    print(den) 

    den= 0   
    print('dcs density') 
    d_l = []        
    for key in range(0,numberOfGraphs):
        den+=comDensity(snapshots[key],dcs)
        d_l.append(comDensity(snapshots[key],dcs))
    # print(d_l)
    # print(min(d_l))
    print(den)  

 # Jaccard values
    print('discovered min jac') 
    
    s1 =set()
    s2 = set()
    
    jac = []
    
    for idx1 in range(0,numberOfGraphs): 
        for idx2 in range(0,numberOfGraphs): 
            if idx1 < idx2 : 
                s1 = set(set_dic[idx1])
    
                s2 = set(set_dic[idx2])
    
            
                val = jaccard_similarity(s1,s2)
                jac.append(val)
    print(min(jac))
    print(np.average(jac))

