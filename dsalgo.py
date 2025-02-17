
"""
The desest subgraph algorithms for single graph and 
The densest common subgraph for multiple graph snapshots.



"""

import networkx as nx
import copy
from datetime import datetime, timedelta
import time
from gurobipy import *



"""
Finds weighted densest subgraph : LP based algorithm
"""

def densest_subgraph_w(G): # assumes G is undirected
    vertices = G.nodes()
    und_edges = G.edges()
    if not und_edges: return 0, []    

    model = Model()

    # Suppress output
    model.params.OutputFlag = 0

    # Add variables
    y = model.addVars(vertices, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="y")
    x = model.addVars(und_edges, lb=0, vtype=GRB.CONTINUOUS, name="x")
    model.update()

    # Size constraint
    model.addConstr(quicksum(y[i] for i in vertices) == 1)

    # Edge constraints
    for v, w in und_edges:
        model.addConstr(x[v, w] <= y[v])
        model.addConstr(x[v, w] <= y[w])

    # Set objective function (average degree)
    # model.setObjective(quicksum(2 * x[v, w]* G.edges[v, w]['weight'] for (v, w) in und_edges))
    model.setObjective(quicksum(x[v, w]* G.edges[v, w]['weight'] for (v, w) in und_edges))
    model.modelSense = GRB.MAXIMIZE

    model.update()
    model.optimize()

    assert model.status == GRB.status.OPTIMAL
    sol = [ v for v in vertices if y[v].x > 0 ]
    
    yis = [y[v].X for v in vertices if y[v].X > 0]
    # print(yis)
    induced = G.subgraph(sol)        
    d =  induced.size(weight="weight") / induced.number_of_nodes()
    # d = sum([induced[u][v]['weight'] for u,v in induced.edges()])  / induced.number_of_nodes()
    # 2.0 * induced.number_of_edges()
    # assert d >= 2.0 * len(und_edges) / len(vertices)
    # print(sol)
    
    obj = model.getObjective()
    # print("Objective = ", obj.getValue())
    return d, induced , sol

"""
IP based algo 

"""


def ip_dcs_sum(alpha, snapshots, comb, eta, A, und_edges, vertices):

    M = pow(10, 5)
    # print("M value = ",M)
    k = len(snapshots)

    model = Model()

    # Suppress output
    model.params.OutputFlag = 0

    # Add variables
    y = model.addVars(vertices, vtype=GRB.BINARY, name="y")
    x = model.addVars(und_edges, vtype=GRB.BINARY, name="x")
    z = model.addVars(und_edges, vtype=GRB.BINARY, name="z")
    model.update()

    # Edge constraints
    for v, w in und_edges:
        # subgraph constraint
        model.addConstr(x[v, w] <= y[v])
        model.addConstr(x[v, w] <= y[w])
        

    model.setObjective(quicksum(1*quicksum(x[v, w] for (v, w) in A[i])
                                for i in range(k)) - (eta*quicksum(y[ii] for ii in vertices)))
    model.modelSense = GRB.MAXIMIZE

    model.update()
    model.optimize()

    # assert model.status == GRB.OPTIMAL
    # if 	GRB.OPTIMAL == 2:
    #     print("Model status = OPTIMAL")
    if GRB.OPTIMAL != 2:
        print("Model status = NOT OPTIMAL")
    obj = model.getObjective()
    # print("Objective = ", obj.getValue())

    yis = [y[v].X for v in vertices if y[v].X > 0]
    subg = [v for v in vertices if y[v].X > 0]
    

    return [obj.getValue(), subg]

"""
Perform binary serach

"""


def ip_based_dcs_sum(alpha1, snapshots, comb):
    eps = 0.01
    # G = snapshots[0]
    n = comb.number_of_nodes()
    up = .5*(n - 1)*len(snapshots)

    vertices = comb.nodes()
    k = len(snapshots)
    low = 0
    # low =  14.0511474609375
    subg = {}
    val = 0
    comb = comb.to_undirected()
    und_edges = list(set([tuple(set(edge)) for edge in comb.edges()]))
    obj_val = 0

    # compute adjacencies
    A = []

    for i in range(k):
        ed_ls = []

        for (a, b) in snapshots[i].edges():
            if (a, b) in und_edges:
               ed_ls.append((a, b))
            elif (b, a) in und_edges:
               ed_ls.append((b, a))
            else:
                print('...')
        A.append(ed_ls)

    itr = 0
    # binary search
    while up - low > eps*low:
        alpha = (up + low) / 2
        # print('alpha {} , low {} , up {}'.format(alpha, low, up))
        
        # obj, solution, s_e, w_e = ip(lam, G, alpha1, w_edges)

        [obj, solution] = ip_dcs_sum(alpha1, snapshots, comb,
                              alpha, A, und_edges, vertices)
        itr += 1
        # print('obj {}  sol {}'.format(obj, solution))
        l = len(solution)

        if l:
            _val = (obj + alpha*l) / l
            # print('alpha: {} obj:  {} len: {}'.format(alpha,_val, l))
        # else:
        #     print('alpha:  {} {}  len = 0'.format(alpha, obj))

        if len(solution) == 0:
            # if not obj:
            up = alpha
        else:
            low = alpha

            if _val > val:
                subg = solution
                val = _val
                obj_val = obj

    print('obj val:  {} subgraph :  {} : itr: {}, 0 {}'.format(val, len(subg), itr,obj_val))

    return [val, subg]

"""
Goldberg's exact max flow algorithm.

Parameters
----------
G: undirected, graph (networkx).

Returns
-------
subg: list, subset of nodes corresponding to densest subgraph.
opt: float, density of subg induced subgraph.

"""
def exact_densest(G: nx.Graph):

    m = G.number_of_edges()
    n = G.number_of_nodes()
    
    low = 0
    up = (n - 1) / 2

    opt = 0
    subg = G.nodes()
    
    if low == up:
        return subg, up
    
    # binary search
    while up - low > 1 / (n *(n -1)): 
        guess = (up + low) / 2
        H = create_flow_network(G, guess)
        
        solution = nx.minimum_cut(H, 's', 't', capacity='capacity')  
        cut = solution[1][0]

        if cut == {'s'}:
            up = guess  
        else:            
            low = guess
            subg = cut
            opt = guess
            
    subg = list(set(subg)&set(G.nodes()))
    
    return subg, opt

def create_flow_network(G, guess):
    m = G.number_of_edges()
    m = int(G.size(weight="weight"))
    G = nx.DiGraph(G)
    
    H = G.copy()
    H.add_node('s')
    H.add_node('t')
    
    for e in G.edges():
        H.add_edge(e[0], e[1], capacity=G.edges[e[0], e[1]]['weight'])

    for v in G.nodes():
        H.add_edge('s', v, capacity=5*m)
        H.add_edge(v, 't', capacity=5*m + 2 * guess - G.in_degree(v, 'weight'))
    return H


def exact_densest_weighted(G: nx.Graph,k):

    m = G.number_of_edges()
    n = G.number_of_nodes()
    
    low = 0
    up = k*(n - 1) 

    opt = 0
    subg = G.nodes()
    
    if low == up:
        return subg, up
    
    # binary search
    while up - low > 1 / (n *(n -1)): 
        guess = (up + low) / 2
        # print('alpha {}'.format(guess))
        H = create_flow_network(G, guess)
        
        solution = nx.minimum_cut(H, 's', 't', capacity='capacity')  
        cut = solution[1][0]

        if cut == {'s'}:
            up = guess  
        else:            
            low = guess
            subg = cut
            opt = guess
            
    subg = list(set(subg)&set(G.nodes()))
    
    return subg, opt

