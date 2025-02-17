

"""
CASE-STUDIES-2
Experiments for DBLP : Smalller Dataset

"""


import numpy as np
import networkx as nx
import copy
import time
import os
import pandas as pd
from dsalgo import *
from find_densest_distinct_sets import *


# read data
filepath = os.path.join("Data","dblp.txt")
header_list =  ['source','target', 'timestamp']
df = pd.read_table(filepath,sep='\t',names=header_list)
# Remove self-loops
df = df[((df['source'] ) != (df['target']))]
keys = df.timestamp.unique()

data_dic = {k: [] for k in keys}
data_dic = {k: [] for k in range(2006, 2016)}


author_cnt = {}
row_indx_lst = []
for _row in df.values:
    if _row[2] < 2016 and _row[2] > 2005:
        
        if _row[0] in author_cnt:
            author_cnt[_row[0]] = author_cnt[_row[0]].union(set({_row[2]}))            
            # print('yes')
        else:
            author_cnt[_row[0]] = set({_row[2]})
            # author_cnt[_row[0]].union({_row[2]})
                
        if _row[1] in author_cnt:
            author_cnt[_row[1]] = author_cnt[_row[1]].union(set({_row[2]}))
        else:
            author_cnt[_row[1]] = set({_row[2]})
            # author_cnt[_row[1]].union({_row[2]})           
            
alst = []
for idx, s in enumerate(author_cnt):   
    print(author_cnt[s])
    if   len(author_cnt[s]) < 4:
        alst.append(s)
        # del author_cnt[idx]

print('size of the removed authors')            
print(len(alst))

data_dic = {k: [] for k in range(2006, 2016)}
df_all= None 
for _row in df.values:
    if _row[2] < 2016 and _row[2] > 2005:
        if (_row[0] not in alst) and  (_row[1] not in alst):
             data_dic[_row[2]].append((_row[0],_row[1]))

snapshots =  []

sub = []
deg = []
nodes = set()

print(len(data_dic.keys()))
#%% snaps
edg_lst = []
avg_nodes = 0
avg_edges = 0
for k, ls in data_dic.items():
    G=nx.Graph()  
    G.add_edges_from(ls, weight = 1)
    snapshots.append(G)
    
    print(len(G.nodes()))
    avg_nodes += len(G.nodes())
    print(len(G.edges()))
    avg_edges += len(G.edges()) 
        
    nodes = nodes.union( set(G.nodes()))
    edg_lst += [e for e in G.edges]
    
numberOfGraphs = len(snapshots)

avg_nodes /= numberOfGraphs 
avg_edges /= numberOfGraphs 

print('avg edges: %s'%avg_edges)

sub = []
deg = []

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


#%%  densest subgraph algo : multi graphs
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

exact_R = exact_densest(G_avg)
# print('subgraph induced by', exact_R[0])
print('density =', exact_R[1])

dcs3 = set(exact_R[0])
print('densest common: unweighted')
print(len(dcs3))
den = 0        
for key in range(0,numberOfGraphs):
    subgrpah_snap = snapshots[key].subgraph(dcs3)
    den+= (subgrpah_snap.number_of_edges()/ len(dcs3) )
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

algo_ver = 0

if algo_ver == 0:
    # hard constraint
    alp = .4
    set_dic =  hard(snapshots, dcs,  alp)  
    print('alp  {}'.format(alp))

elif algo_ver == 1: 
    # soft contraint - 1
    lam = .8
    [s1,set_dic1] = soft_1_1(snapshots, ld, lam)
    [s2, set_dic2] = soft_1(snapshots, dcs, lam)
    print('lam  {} s1: {} s2: {}'.format(lam,s1,s2))
    
    if s1 > s2 :
        set_dic = set_dic1
    else:
        set_dic = set_dic2
 
elif algo_ver == 2: 
    # soft contraint - 2
    lam = 5
    set_dic = soft_2(snapshots, dcs, lam, nodes)
    print('lam  {}'.format(lam))

else:
    print('no algo found')
    

print('sets density') 
den= 0     
d_l = []      
# read data
filepath = os.path.join("Data","dblp_Authors.txt")
header_list =  ['id','name']
des = pd.read_table(filepath,sep='\t',names=header_list)
# des.loc[des['id'] == 0].values[0][1]
print('sets density') 
den= 0     
d_l = []      
for key, val in set_dic.items(): 
    print(key+2006)
    print(len(val))
    s1= val
    inter = dcs.intersection(s1)
    extra_add = s1- dcs
    print()
    print('added')
    for i in extra_add:
        str1 = des.loc[des['id'] == int(i)].values[0][1]
        str2 = str1.split(" ")
        str2[0] = [*str2[0]][0] + '.'
        new_st = ""
        # a = 0
        for j in str2:
            # a+=1
            new_st= new_st + j
            # if a < (len(str2) - 1):
            # new_st = new_st + " "
            
        # print(des.loc[des['id'] == int(i)].values[0][1], end =", ")
        print(new_st, end =", ")
    
    removed = dcs -s1
    print()
    print('removed')
    for i in removed:
        str1 =des.loc[des['id'] == int(i)].values[0][1]
        str2 = str1.split(" ")
        str2[0] = [*str2[0]][0] + '.'
        new_st = ""
        # a=0
        for j in str2:
            # a+=1
            new_st= new_st + j
            # if a < (len(str2) - 1):
            # new_st = new_st + " "
            
        # print(des.loc[des['id'] == int(i)].values[0][1], end =", ")
        print(new_st, end =", ")
        # print(des.loc[des['id'] == int(i)].values[0][1], end =", ")
        
    den+=comDensity(snapshots[key],val)
    d_l.append(comDensity(snapshots[key],val))
    print('--------------')
print(d_l)
print(min(d_l))
print(den)  

den= 0   
print('dcs density') 
d_l = []        
for key in range(0,numberOfGraphs):
    
    den+=comDensity(snapshots[key],dcs)
    d_l.append(comDensity(snapshots[key],dcs))
print(d_l)
print(min(d_l))
print(den) 

print('size {} set {}'.format(len(dcs),dcs))
for i in dcs:
    # print(des.loc[des['id'] == int(i)].values[0][1], end =", ")
    str1 = des.loc[des['id'] == int(i)].values[0][1]
    str2 = str1.split(" ")
    str2[0] = [*str2[0]][0] + '.'
    new_st = ""
    # a=0
    for j in str2:
        # a+=1
        new_st= new_st + j
        # if a < (len(str2) - 1):
        # new_st = new_st + " "
        
    # print(des.loc[des['id'] == int(i)].values[0][1], end =", ")
    print(new_st, end =", ")
        

 #%% Jaccard values

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
            # print(val)
            jac.append(val)
print(min(jac))



