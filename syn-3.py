
"""
Synthetic-dataset-3
"""


import random
import networkx as nx
import pandas as pd
import os
from dsalgo import *
from find_densest_distinct_sets import *
import pickle

# dataset params

seed = 30528443483785
numberOfGraphs =5

denseGProb = .06
sparseGProb =.005
crossGProb = .002

denseVertices = 350
sparseVertices = 3500

noise_prob = [.01*random.randint(0, 9) for i in range (0, numberOfGraphs)]

def generate_data(numberOfGraphs,denseGProb,sparseGProb,crossGProb,noise_prob, denseVertices,sparseVertices,seed):
    edg_lst = []
    nodes = set()
    snapshots =  []
    
    DenseGraphSets = [list(range(denseVertices)) for i in range (0, numberOfGraphs)]
    SparseGraphSets = [list(range(denseVertices, sparseVertices)) for i in range (0, numberOfGraphs)]
    
    df_all = None
    column_names = ["source", "target","snap" ]
    df_all = pd.DataFrame(columns=column_names)
        
    for k in range (0, numberOfGraphs):    
        print('snapshot : {}'.format( k))
        sparse_set = range(denseVertices,sparseVertices)
        length =  int(len(sparse_set)*noise_prob[k])
        
        noisy_set = random.sample(sparse_set,length)
        DenseGraphSets[k] += noisy_set
        
        SparseGraphSets[k] = list(set(SparseGraphSets[k]) -  set(noisy_set))
        
        U=nx.Graph()
        
        d_edges= 0
        # dense edges
        for i in DenseGraphSets[k]:
            for j in DenseGraphSets[k]:
                if (i < j):
                    R = random.random()
                    if (R < denseGProb):
                        U.add_edge(i, j)
                        d_edges += 1
                        df = pd.DataFrame([[i, j,k]], columns=column_names)            
                        df_all = pd.concat([df_all, df], ignore_index=True)
        print('dense edges: {}'.format( d_edges))
        
        s_edges= 0
        # sparse edges
        for i in SparseGraphSets[k]:
            for j in SparseGraphSets[k]:
                if (i < j):
                    R = random.random()
                    if (R < sparseGProb):
                        U.add_edge(i, j)
                        s_edges += 1
                        df = pd.DataFrame([[i, j,k]], columns=column_names)            
                        df_all = pd.concat([df_all, df], ignore_index=True)
        print('sparse edges: {}'.format( s_edges))
        
        c_edges = 0
        # cross edges
        for i in DenseGraphSets[k]:
            for j in SparseGraphSets[k]:
                if (i < j):
                    R = random.random()
                    if (R < crossGProb):
                        U.add_edge(i, j)
                        c_edges+= 1
                        df = pd.DataFrame([[i, j,k]], columns=column_names)            
                        df_all = pd.concat([df_all, df], ignore_index=True)
        print('cross edges: {}'.format( c_edges))

        edg_lst += [e for e in U.edges]
        nodes = nodes.union( set(U.nodes()))

        print(len(U.nodes()))
        print(len(U.edges()))
        print(length)
        
        snapshots.append(U)
            
    df_all = df_all[['source', 'target','snap']]
    
    # Save as .csv file
    df_all.to_csv('./Data/syn-sen-5.csv')
    # df_all.to_csv('./Data/t2.txt', sep=' ', index=False, header=None)
    
    with open('./Data/d-sen-5', 'wb') as fp:
        pickle.dump(DenseGraphSets, fp)
    with open('./Data/s-sen-5', 'wb') as fp:
        pickle.dump(SparseGraphSets, fp)
        
    return [snapshots, edg_lst, nodes,DenseGraphSets,SparseGraphSets ]

# [snapshots, edg_lst, nodes,DenseGraphSets,SparseGraphSets ] = generate_data(numberOfGraphs,denseGProb,sparseGProb,crossGProb,noise_prob, denseVertices,sparseVertices,seed)


#%% 

""" read data """

DenseGraphSets = []
with open ('./Data/d-sen-5', 'rb') as fp:
    DenseGraphSets = pickle.load(fp)
    
SparseGraphSets= []

with open ('./Data/s-sen-5', 'rb') as fp:
    SparseGraphSets = pickle.load(fp)

# read data
filepath = os.path.join("Data", "syn-sen-5.csv")
df = pd.read_csv(filepath)

keys = df.snap.unique()
print(keys)

data_dic = {k: [] for k in keys}

for _row in df.values:
        data_dic[_row[3]].append((_row[1],_row[2]))

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
 #%%
print('dcs density LP') 
d, induced , dcs1 = densest_subgraph_w(G_avg)
print(len(dcs1))
den = 0        
for key in range(0,numberOfGraphs):
    subgrpah_snap = snapshots[key].subgraph(dcs1)
    den+= (subgrpah_snap.number_of_edges()/ len(dcs1) )
print(den)
dcs = set(dcs1)
#%%  DSC and d_ind correct
# comb = snapshots[0]
# for i in range(1, len(snapshots)):

#     H = snapshots[i]
#     comb = nx.compose(comb, H)

# nodes = comb.nodes()

# [obj1, subg1] =ip_based_dcs_sum(0, snapshots, comb)

# densities = []
# for i in range(len(snapshots)):
#     G = snapshots[i]
#     sub_g = G.subgraph(subg1)
#     den = sub_g.number_of_edges()/len(subg1)
#     densities.append(den)
#     # print(den, end=", ")
# print('IP-based')
# print(sum(densities))

# dcs = set(subg1)
  #%%
ld = []
local_sum_den = 0
for snap in snapshots:
    d, induced ,  sol = densest_subgraph_w(snap)
    local_sum_den += d
    ld.append(set(sol))
print('local sum denisty')    
print(local_sum_den )  

  #%% Algorithms
algo_ver = 1

if algo_ver == 0:
    # hard constraint
    alp = .7
    set_dic =  hard(snapshots, dcs,  alp) 
    print('alpha  {}'.format(alp))

elif algo_ver == 1: 
    # soft contraint - 1
    lam = .7

    [s1,set_dic1] = soft_1_1(snapshots, ld, lam)
    [s2, set_dic2] = soft_1(snapshots, dcs, lam)
    print('lam  {} s1: {} s2: {}'.format(lam,s1,s2))
    
    if s1 > s2 :
        set_dic = set_dic1
    else:
        set_dic = set_dic2
    print('lam  {}'.format(lam))
elif algo_ver == 2: 
    # soft contraint - 2
    lam = 0.01
    set_dic = soft_2(snapshots, dcs, lam, nodes)
    print('lam  {}'.format(lam))

else:
    print('no algo found')
    
 #     Results

print('sets density') 
den= 0     
d_l = []      
for key, val in set_dic.items():   
    # print('size {} set {}'.format(len(val),val))
    print('size {} '.format(len(val)))
    den+=comDensity(snapshots[key],val)
    d_l.append(comDensity(snapshots[key],val))
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


den= 0    
print('ground truth density') 
d_l = []        
for key in range(0,numberOfGraphs):
    den+=comDensity(snapshots[key],set(DenseGraphSets[key]))
    d_l.append(comDensity(snapshots[key],set(DenseGraphSets[key])))
print(d_l)
print(min(d_l))
print(den)

 

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
print('discovered avg jac') 
print(np.average(jac))

print('ground truth min jac') 

s1 =set()
s2 = set()

jac = []

for idx1 in range(0,numberOfGraphs): 
    for idx2 in range(0,numberOfGraphs): 
        if idx1 < idx2 : 
            s1 =  set(DenseGraphSets[idx1])

            s2 = set(DenseGraphSets[idx2])

        
            val = jaccard_similarity(s1,s2)
            jac.append(val)
print(min(jac))


val = 0

for idx in range(0,numberOfGraphs):      
    s1 = set(set_dic[idx])

    s2 = set(DenseGraphSets[idx])


    val += jaccard_similarity(s1,s2)
print('average jacaard')
print(val/numberOfGraphs)  
