

"""
Algorithms (Naive implementation without priority queues)

Hard
Soft-I
Soft-II

"""

import copy
import time
import math
import heapq as heap
import numpy as np
from operator import attrgetter
    
""" Vertex structure"""
class Vertex:
    def __init__(self, data):
            self.key = data
            self.int= []
            self.un = []
            self.deg = 0   
            self.den = 0
            self.cons = 0
            self.tot_gain = 0
            self.ind = []
            self.pattern = []
            
    def setTotGain(self):
        self.tot_gain =  self.den + self.cons
        
    def setInt(self,i):
        self.int[i] -= 1
        
    def setUn(self,i):
        self.un[i] -= 1
        
    def setDen(self,diff):
        self.den += diff 
        
    def setDeg(self,deg):
        self.deg = deg
        
""" Snapshot structure
Includes pattern clusters and sum contraint value for each pattern cluster
"""
class Snapshot:
   def __init__(self, data):
       self.key = data
       self.m = 0
       self.n = 0
       # Store the patten clusters in a dictionary 
        # { 0010 : [ 1, 2 ], 1011 : [ 3, 4 ] }
       self.pat_dic = {} 
       #store sum constraint value specific to each pattern
       self.pat_val = {}
       # #store priority queues for each pattern
       # self.pat_pq = {}
       
       
   def AddPattern (self, pattern, node ) :
       if pattern in self.pat_dic  :
           self.pat_dic[pattern].append(node)
           self.pat_val[pattern] = 0
       else:
           self.pat_dic[pattern] = [node]
           if pattern not in self.pat_val :
                self.pat_val[pattern] = 0
           
   def removeElement(self,pattern,node):
       self.pat_dic[pattern].remove(node)
       
   def deletePattern(self,pattern):
       del self.pat_dic[pattern]
       del self.pat_val[pattern]
           
   def updatePattern(self,  node, snap, max_snap,pattern) :
       
           lst = self.pat_dic[pattern]
           # currnt_pat_val = self.pat_val[pattern]
           
           if node in  lst:
                new_pattern= list(pattern)
                lst.remove(node)
                
                if len(lst) == 0:
                    self.pat_dic.pop(pattern)
                    self.pat_val.pop(pattern)
                new_pattern[max_snap] = 0
                
                pattern = tuple(new_pattern)
                if pattern in self.pat_dic  :
                    self.pat_dic[pattern].append(node)
                    self.pat_val[pattern] = 0
                else:
                    self.pat_dic[pattern] = [node]
                    if pattern not in self.pat_val :
                        self.pat_val[pattern] = 0
                       
           return new_pattern       
                                    
   def Display_Patterns(self) :
        for item in self.pat_dic.items() :
            print(item)
            
   def Display_Patter_val(self) :
        for item in self.pat_val.items() :
            print(item)
       
   def setM(self, m):
       self.m = m
   
   def setN(self, n):
       self.n = n
    
"""
Computes density of a subgraph of given vertex set
"""


def comDensity(G, S):
    val = 0
    if len(S) != 0:
        val = len(G.subgraph(S).edges) / len(S)
        return val
    return 0


"""
Computes pairwise Jaccard similarity
"""


def jaccard_similarity(A, B):
    similarity = 0

    nominator = A.intersection(B)
    denominator = A.union(B)

    if len(denominator) != 0:
        similarity = len(nominator)/len(denominator)

    return similarity


"""
Computes summation of pairwise Jaccard similarities 
"""


def sum_con(set_dic):
    numberOfGraphs = len(set_dic)
    _sum = 0
    # print(set_dic[0])
    for idx1 in range(numberOfGraphs):
        for idx2 in range(numberOfGraphs):
            if idx1 < idx2:

                s1 = set(set_dic[idx1])
                s2 = set(set_dic[idx2])

                val = jaccard_similarity(s1, s2)
                _sum += val

    return _sum

"""
Check Jac. Constraint efficient
"""
def checkJacCon(I,U, alp, ind, v, set_dic, isRem):
    isTrue = True
    s1 = set()
    
    
    _numberOfGraphs = len(set_dic)
    for i in range(_numberOfGraphs):
        if ind != i:
            s1 = set_dic[i]
            val = 1
            
            if isRem: 
                if v in s1:
                    if U[i, ind] != 0:
                        val = (I[i, ind] - 1) / U[i, ind]
    
                else:
                    if U[i, ind] > 1:
                        val = (I[i, ind]) / (U[i, ind] - 1)
            else:  
                if v in s1: 
                    if U[i, ind] != 0:
                        val = (I[i, ind] + 1) / U[i, ind]
                
                else:
                    val = (I[i, ind]) / (U[i, ind] + 1)
                    
            if val < alp:
                
                
                isTrue = False
                return isTrue
                break
            

    return isTrue

"""
Computes the difference in pairwise Jaccard values  efficiently
if vertex v is removed from snapshot ind
"""


def sum_con_efficient(I, U, set_dic, ind, v ):
    _numberOfGraphs = len(set_dic)
    _sum = 0
    
    for i in range(_numberOfGraphs):

        if i != ind:

            s1 = set_dic[i]
            
            init = I[i, ind]/U[i, ind]
            val1 = ((I[i, ind] - 1) / U[i, ind]) - init
            val2 = (I[i, ind] / (U[i, ind] - 1)) - init
            
            
            if v in s1:
                    _sum += val1

            else:
                if U[i, ind] != 1:
                    _sum += val2
    # print(_sum)                
    return _sum

"""
Computes  SUM CONSTRAINTs based on vertex pattern
"""
def sum_con_pattern(I, U, ind, pattern ):
    _numberOfGraphs = len(pattern)
    _sum = 0
    
    for i in range(_numberOfGraphs):

        if i != ind:

            # s1 = set_dic[i]
            
            init = I[i, ind]/U[i, ind]
            val1 = ((I[i, ind] - 1) / U[i, ind]) - init
            val2 = (I[i, ind] / (U[i, ind] - 1)) - init
            
            
            # print(init)
            
            # if v in s1:
            if pattern[i]:
            # if  ((i,v) in bool_dic) :
            # if bool_dic[i,v]:
                    _sum += val1

            else:
                if U[i, ind] != 1:
                    _sum += val2
    # print(_sum)                
    return _sum

"""
Update I & U
"""

def updateIU(I, U, set_dic, ind, v, isRem =  True):
    _numberOfGraphs = len(set_dic)
    # _sum = 0
    # print(I)
    # print(U)

    for i in range(_numberOfGraphs):

        if i != ind:

            s1 = set_dic[i]
           #  un = len(s1.union(S))
            # print(type(s1))
            if isRem:
                if v in s1:
                    I[i, ind] = (I[i, ind] - 1)
                    I[ind, i] = I[i, ind]
    
                else:
                    U[i, ind] = (U[i, ind] - 1)
                    U[ind, i] = U[i, ind]
            else:
                if v in s1:
                    I[i, ind] = (I[i, ind] + 1)
                    I[ind, i] = I[i, ind]

                else:
                    U[i, ind] = (U[i, ind] + 1)
                    U[ind, i] = U[i, ind]

"""
Computes  f = q() + lam()
"""


def com_q(snapshots, set_dic, lam):
    # print('sets density')
    den = 0
    d_l = []
    for key, val in set_dic.items():
        # print('size {} set {}'.format(len(val),val))
        # print('size {} '.format(len(val)))
        den += comDensity(snapshots[key], val)
        d_l.append(comDensity(snapshots[key], val))
    # print(d_l)
    # print(min(d_l))

    den += (lam*sum_con(set_dic))

    return den




"""
Check Jac. Constraint
"""


def checkJacConstraint(alp, snapshots, i, S, set_dic):
    isTrue = True
    s1 = set()
    s2 = set()

    for idx1, snap_G1 in enumerate(snapshots):
        for idx2, snap_G2 in enumerate(snapshots):
            if idx1 < idx2:

                if idx1 != i:
                    s1 = set(set_dic[idx1])
                else:
                    s1 = S

                if idx2 != i:
                    s2 = set(set_dic[idx2])
                else:
                    s2 = S

                val = jaccard_similarity(s1, s2)
                # print(val)

                if val < alp:
                  isTrue = False
                  return isTrue
                  # break

    # print(isTrue)
    return isTrue


"""
Check whether there is any non-empty sets remaining

"""


def checkEmpty(set_itr):
    isEmpty = True
    for k, itm in set_itr.items():
        if len(itm) > 0:
            # print(type(itm))
            isEmpty = False
            return isEmpty
    return isEmpty


"""
Hard constraint algo - ALGO 0

"""


def hard(snapshots, dcs,  alp):
    list_of_sets = [dcs for itm in snapshots]
    set_dic = {}
    numberOfGraphs = len(snapshots)
    keys = range(numberOfGraphs)
    values = list_of_sets
    set_dic = dict(zip(keys, values))

    _prev_val = pow(10, 3)
    _curr_val = 0
    _itr = 0

    start_time = time.time()
    
    #initial intersections and unions between all pairs of sets
    I = {}
    U = {}

    for i in range(numberOfGraphs):
        for j in range(numberOfGraphs):
            if i < j:
                val = len(dcs)
                I[i, j] = val
                I[j, i] = val
                U[i, j] = val
                U[j, i] = val

    while round(_prev_val, 1) != round(_curr_val, 1):
        print("iteration no........... %d " % (_itr+1))
        _prev_val = _curr_val

        for idx, snap_G in enumerate(snapshots):
            S_i = set(snap_G.nodes())

            for v in S_i:
                d1 = 0
                d2 = 0
                
                S = copy.deepcopy(set_dic[idx])

                current_g = snap_G.subgraph(S)
                m = len(current_g.edges())
                n = len(current_g.nodes())

                if S and (len(S) > 0) and (v in S):
                    # S.remove(v)
                    
                    if checkJacCon(I,U, alp, idx, v, set_dic, True):
                        deg = current_g.degree(v)
                        
                        if n != 1:
                            d1 = -( (m - deg) / (n - 1)) + (m / n)

                # S = copy.deepcopy(set_dic[idx])

                if S and (v not in S):
                    S.add(v)

                    if checkJacCon(I,U, alp, idx, v, set_dic, False):
                        deg = (snap_G.subgraph(S.union({v}))).degree(v)
                        
                        if n != 0:
                            d2 = - ((m + deg) / (n + 1)) + (m / n)

                d = min(d1, d2)

                if d < 0:
                    isR = d1 < d2

                    if isR:
                        set_dic[idx] = set_dic[idx] - {v}
                        # set_dic[idx].remove(v)
                        updateIU(I, U, set_dic, idx, v, True )
                    else:
                        set_dic[idx] = set_dic[idx].union({v})
                        # set_dic[idx].add(v)
                        updateIU(I, U, set_dic, idx, v, False)

                    _curr_val -= d

        print('current val: {} - prev: {}'.format(_curr_val, _prev_val))
        _itr += 1

    print("--- %s seconds ---" % (time.time() - start_time))
    return set_dic

# intialize snapshots for iterative algorithm when starting off with the gfull set of the vertices
def _initialize(snapshots_org,set_dic,lam,I,U,idx):

    # dictionary to store vertex objects
    v_count_dic = {}

    
     # initial snapshots :copy original snapshots
    snapshots = copy.deepcopy(snapshots_org)
        
    # vertex counters initialization
    # for idx, snap in enumerate(snapshots):
    snap = snapshots[idx]
            
    snapObj =  Snapshot(idx)    
    m = len(snap.edges())
    n = len(snap.nodes())
    snapObj.m = m
    snapObj.n = n
    
    for i in snap:

        obj = Vertex(i)
        
        # set vertex counters
        obj.deg = snap.degree(i)
        
        _pattern = []

        for ii, (k, lst1) in enumerate(set_dic.items()):
                if ii == idx:
                    _pattern.append(-1)
                elif i in lst1:
                    _pattern.append( 1)
                else:
                    _pattern.append( 0)
        obj.pattern = _pattern 
        snapObj.AddPattern(tuple(_pattern), i)
        
        v_count_dic[i] = obj
       

    # S_i = set_dic[idx]
    S_i = snap
    
    m = snapObj.m
    n = snapObj.n

    for v in S_i:
        
        deg = v_count_dic[v].deg
        if n != 1:
            den = (m - deg) / (n - 1) - (m / n) 
        else:
            den = 0
        # if v in   v_count_dic: 
        v_count_dic[v].den  = den
                            

    pattern_dic = snapObj.pat_val
    
    for _pattern, val in pattern_dic.items():
        # print(_pattern)
        pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
            # print(pattern_dic[_pattern])
    return v_count_dic, snapObj      




"""
Soft constraint algo - 1 : density based,  peeling one snapshot fully

"""


def soft_1(snapshots, dcs, lam):
    numberOfGraphs = len(snapshots)
    list_of_sets = [dcs for itm in snapshots]
    set_dic = {}
    keys = range(numberOfGraphs)
    values = list_of_sets
    set_dic = dict(zip(keys, values))

    d = com_q(snapshots, set_dic, lam)
    print('current density: {}'.format(d))
 
    start_time = time.time()

    #initial intersections and unions between all pairs of sets
    I = {}
    U = {}
    
    for i in range(numberOfGraphs):
        for j in range(numberOfGraphs):
            if i < j:
            
                I[i, j] = len(set_dic[i].intersection(set_dic[j]))
                U[i, j] = len(set_dic[i].union(set_dic[j]))               
                I[j, i] = I[i, j]
                U[j, i] = U[i, j]  
    
    _prev_val =  [-1 for i in range(numberOfGraphs)]
    _curr_val =  [0 for i in range(numberOfGraphs)]
    
    isImp = True
    _itr = 0
    
    # iterative algo
    while isImp :
        print("iteration no........... %d " % (_itr+1))
        _prev_val = copy.deepcopy(_curr_val)
        isImp = False
        
        for idx, snap_G in enumerate(snapshots):
                    
            # initialize with  set of all nodes 
            S_i = copy.deepcopy(set(snap_G.nodes()))
            
            # current set
            current_set = set_dic[idx]
            
            _max  = 0 # max change
            _vals = 0 # current change
            
            # update I and U 
            for i in range(numberOfGraphs):
                if i != idx:                        
                    I[i, idx] = len(
                        set(snap_G.nodes()).intersection(set_dic[i]))
                    U[i, idx] = len(
                        set(snap_G.nodes()).union(set_dic[i]))

                    I[idx, i] = I[i, idx]
                    U[idx, i] = U[i, idx]     
                    
            v_count_dic, snap_obj = _initialize(snapshots,set_dic,lam,I,U,idx)  
            snapshots1 = copy.deepcopy(snapshots[idx])
            ################################ SPEED-UP############
            while len(S_i) > 0:
                _maxx = -math.inf
                max_snap = idx
                max_v = -1

                _snap = idx                                
                patt_dic =  snap_obj.pat_dic
                patt_val = snap_obj.pat_val
    
                for pat in patt_val.keys():  
                    _max_pat = -math.inf
                    _max_pat_v= None
                    
                    #search for the vertex with max density component within a pattern
                    for v in patt_dic[pat]:
                        gain_v = v_count_dic[v].den
                        
                        if gain_v > _max_pat:
                            _max_pat = gain_v
                            _max_pat_v = v
                    
                    # tot_gain of the max density candidate chosen from each pattern
                    tot_gain = v_count_dic[_max_pat_v].den  + patt_val[pat]         
                        
                    if tot_gain > _maxx:
                        _maxx = tot_gain
                        
                        max_snap, max_v = _snap,_max_pat_v 
                        max_pattern = pat
                            
                # print('val: {} {}'.format(max_snap, max_v))        
                    
                # current snapshots
                chg_snap = snapshots1
                # find neighbours of node v
                neighbour_nodes  = [n for n in chg_snap.neighbors(max_v)]
                
                diff =    _maxx

                #update SNAPSHOT                
                snapshots1.remove_node(max_v)            
                
                # change upon a vertex removal
                _vals += diff  

                # # # update intersection and union matrices
                updateIU(I, U, set_dic, max_snap, max_v)
                
                # update m and n counters after max_v is removed from max_snap
                snap_obj.m = snap_obj.m - v_count_dic[max_v].deg
                snap_obj.n = snap_obj.n -  1
                        
                #remove max_v from max_snap 
                max_pat_dic_lst =  snap_obj.pat_dic[max_pattern]
                max_pat_dic_lst.remove(max_v)
                
                if len(max_pat_dic_lst) == 0:
                    snap_obj.pat_dic.pop(max_pattern)
                    snap_obj.pat_val.pop(max_pattern)
                                
                # remove max_v
                del v_count_dic[max_v]
                
                # remove max_v        
                S_i.remove(max_v)
                
                  # maintain GLOABL diff
                if _vals  > _max:
                    current_set = copy.deepcopy(S_i) 
                    _max = _vals
                        
                m = snap_obj.m
                n = snap_obj.n
                
                for v in S_i:
                    obj  = v_count_dic[v]  
        
                    # update degree of neigbour nodes of v
                    if v in neighbour_nodes:
                            # update degree
                            obj.deg -= 1                                
                    
                    deg = obj.deg
                    
                    # if deg != snapshots1.degree[v]:
                    #     print('d {} {}'.format(deg, snapshots1.degree[v]))
                    # update density
                    if n != 1:
                        obj.den  = ((m - deg) / (n - 1)) - (m / n)  
                    else:
                        obj.den = 0                                          
        
                #update pattern based sum contraint gains: since I & U are modified now

                pattern_dic = snap_obj.pat_val
    
                for _pattern, val in pattern_dic.items():
                    # print(_pattern)
                    pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
                    
            ###############################NAIVE########
            # start peeling process greedily
            
            # current_g = copy.deepcopy(snap_G)
            # while len(S_i) > 0:
            #     pq = []

            #     m = len(current_g.edges())
            #     n = len(current_g.nodes())

            #     # iterate over all remaining vertices
            #     for v in S_i:
                    
            #         if n != 1:
            #             deg = current_g.degree(v)
            #             sc = (m - deg) / (n - 1) - (m / n)
            #         else:
            #             sc = 0
            #         sc += (lam * sum_con_efficient(I, U, set_dic,
            #                                         idx, v))
                    
            #         heap.heappush(pq, (-sc, v))
                
            #     # pick one vertex to peel
            #     _, argmin = heap.heappop(pq)
            #     #  change upon a vertex removal
            #     _vals += (-_)
                
                
                                
            #     # peel one vertex
            #     S_i.remove(argmin)
            #     current_g.remove_node(argmin)
                
            #     updateIU(I, U, set_dic, idx, argmin)
                
            #     # identify which snapshot produced the best result during peeling process
            #     if _vals > _max:
            #         _max = _vals
            #         current_set = copy.deepcopy(S_i)
            
            ################################
            
            # if the resultant of the peeling, gained an improvement         
            if _curr_val[idx] <  _max:
                _curr_val[idx] = _max
                
            
            if  _curr_val[idx] >  _prev_val[idx]:     
                isImp = True
                set_dic[idx] = current_set
                
   
                         
                
        # print('current val: {} - prev: {}'.format(_curr_val, _prev_val))
        _itr += 1

    print("--- %s seconds ---" % (time.time() - start_time))
    print(' %s iterations'% (_itr))
    print(com_q(snapshots, set_dic, lam))
    return set_dic




"""
Soft constraint algo - 2 : density based,  over all the snapshots

"""


def soft_2(snapshots, dcs, lam, nodes):
    list_of_sets = [set(copy.deepcopy(itm.nodes())) for itm in snapshots]
    numberOfGraphs = len(snapshots)
    set_dic = {}
    keys = range(numberOfGraphs)
    values = list_of_sets
    set_dic = dict(zip(keys, values))

    # set iterations
    list_of_sets_itr = [set(copy.deepcopy(itm.nodes())) for itm in snapshots]
    set_itr = {}
    keys = range(numberOfGraphs)
    values = list_of_sets_itr
    set_itr = dict(zip(keys, values))

    d = com_q(snapshots, set_dic, lam)
    print('current density: {}'.format(d))
    
    start_time = time.time()  
    
    # dictionary to store vertex objects
    v_count_dic = {}
    #dictionary to store snap objects
    snap_dic = {}
    
    # vertex counters initialization
    for idx, snap in enumerate(snapshots):

        I = []
        U = []
        
        for j in range(numberOfGraphs):
            if idx != j :
                I.append( len(set_dic[idx].intersection(set_dic[j])))
                # print('{} {}'.format(len(set_dic[idx].intersection(set_dic[j])), len(set_dic[idx].union(set_dic[j]))))
                U.append(len(set_dic[idx].union(set_dic[j])))
            else:
                I.append( len(set_dic[idx]))
                U.append( len(set_dic[idx]))
            
        snapObj =  Snapshot(idx)    
        m = len(snap.edges())
        n = len(snap.nodes())
        snapObj.m = m
        snapObj.n = n
        
        for i in snap:
            
            obj = Vertex((idx,i))
            
            # set vertex counters
            obj.deg = snap.degree(i)
            
            _pattern = []

            for ii, (k, lst1) in enumerate(set_dic.items()):
                if ii == idx:
                    _pattern.append(-1)
                elif i in lst1:
                    _pattern.append( 1)
                else:
                    _pattern.append( 0)
            obj.pattern = _pattern 
            # snapObj.AddPattern(tuple(_pattern), obj.key)
            snapObj.AddPattern(tuple(_pattern), i)
            
            v_count_dic[idx,i] = obj
        snap_dic[idx] = snapObj
        # snapObj.Display_Patterns()     
        # snapObj.Display_Patter_val()
    # print(v_count_dic[1,110])
         
    delta = 0
    maxx = 0

    # initial snapshots
    snapshots1 = copy.deepcopy(snapshots)
    # initialize I and U matrices
    I = {}
    U = {}
    for i in range(numberOfGraphs):
        for j in range(numberOfGraphs):
            if i < j :
                I[i, j] = len(set_dic[i].intersection(set_dic[j]))
                I[j, i] = I[i, j]  
                U[i, j] = len(set_dic[i].union(set_dic[j]))
                U[j, i] = U[i, j] 
                
    sum_cons_lst = []
    
    for idx1 in range(numberOfGraphs):
        _sum = 0
        for idx2 in range(numberOfGraphs):
            if idx1 != idx2:

                _sum += (lam*I[idx1,idx2]/U[idx1,idx2])
        sum_cons_lst.append(_sum  )  
        
    # initialization         
    for idx, snap_G in enumerate(snapshots1):

        if len(set_itr[idx]) > 0:

            S_i = set_itr[idx]
            
            snap = snap_dic[idx]
            
            m = snap.m
            n = snap.n

            for v in S_i:
                
                deg = v_count_dic[idx,v].deg
                
                if n != 1:
                    den = (m - deg) / (n - 1) - (m / n) 
                else:
                    den = 0

                v_count_dic[idx,v].den  = den
                            
    # for item in snap_dic.items():
    for idx in range(numberOfGraphs):
        snapObj = snap_dic[idx]
        pattern_dic = snapObj.pat_val
        
        for _pattern, val in pattern_dic.items():
            # print(_pattern)
            pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
            # print(pattern_dic[_pattern])
                   
    # iterate through all the nodes
    while not checkEmpty(set_itr):
        # find vertex with max total gain   (GREEDY DIFF (local))              
        _max = -math.inf
        max_snap = -1
        max_v = -1
        # print(numberOfGraphs)
        for _snap in range(numberOfGraphs):
                        
            snapObj = snap_dic[_snap] 
            patt_dic =  snapObj.pat_dic
            patt_val = snapObj.pat_val

            for pat in patt_val.keys():  
                    _max_pat = -math.inf
                    _max_pat_v= None
                    
                    #search for the vertex with max density component within a pattern
                    for v in patt_dic[pat]:
                        gain_v = v_count_dic[_snap,v].den
                        
                        if gain_v > _max_pat:
                            _max_pat = gain_v
                            _max_pat_v = v

                    tot_gain = v_count_dic[_snap,_max_pat_v].den  + patt_val[pat]
                        
                    if tot_gain > _max:
                        _max = tot_gain
                        
                        max_snap, max_v = _snap,_max_pat_v 
                        max_pattern = pat
                    
        # print('val: {} {}'.format(max_snap, max_v))        
            
        # current snapshots
        chg_snap = snapshots1[max_snap]
         # find neighbours of node v
        neighbour_nodes  = [n for n in chg_snap.neighbors(max_v)]
        
        diff =    _max
        # remove the vertex with the highest score
        #update SETS in each iteration
        set_itr[max_snap] = set_itr[max_snap] - {max_v}
        #update SNAPSHOT
        
        snapshots1[max_snap].remove_node(max_v)            
        
        delta += diff  
        # print('{} {}'.format(diff, obj.tot_gain))
        new_sets = copy.deepcopy(set_itr)
        
        # maintain GLOABL diff
        if delta  > maxx:
            set_dic = new_sets
            maxx = delta
       
        # U_old = copy.deepcopy(U)
        # I_old = copy.deepcopy(I)
        
        # # # update intersection and union matrices
        updateIU(I, U, set_dic, max_snap, max_v)
        
        # new_pat_lst = []
        
        #Update patterns of max_v from other snapshots
        for idx, snap_G in enumerate(snapshots1):
            
            if idx != max_snap:

                if max_v in snap_G:                    

                    snap_obj = snap_dic[idx]

                    pattern = tuple(v_count_dic[idx,max_v].pattern)
                    
                    # i = max_snap
                    # _sum = (-lam*(( - 1)/U_old[i, idx]  )) + (lam*((I[i, idx]/(U[i, idx]- 1)) - (I[i, idx]/U[i, idx])))
                    
                    new_pattern= snap_obj.updatePattern(max_v,idx,max_snap,pattern) 
                    v_count_dic[idx,max_v].pattern = new_pattern
                #     new_pat_lst.append(new_pattern)
                # else:
                #     new_pat_lst.append([0])
            # else:
            #      new_pat_lst.append([0])
                 
        #  max_snap snap Obj                           
        maxsnapObj = snap_dic[max_snap] 
        
        # update m and n counters after max_v is removed from max_snap
        maxsnapObj.m = maxsnapObj.m - v_count_dic[max_snap,max_v].deg
        maxsnapObj.n = maxsnapObj.n -  1
                
        #remove max_v from max_snap 
        max_pat_dic_lst =  maxsnapObj.pat_dic[max_pattern]
        max_pat_dic_lst.remove(max_v)
        
        if len(max_pat_dic_lst) == 0:
            maxsnapObj.pat_dic.pop(max_pattern)
            maxsnapObj.pat_val.pop(max_pattern)
        
        #     print('tot {} {}'.format(obj.tot_gain, diff))
        # remove max_v
        del v_count_dic[max_snap,max_v]
        
        # update new density gains in same snapshot        
        S_i = set_itr[max_snap]
                
        m = maxsnapObj.m
        n = maxsnapObj.n
        
        for v in S_i:
            obj  = v_count_dic[max_snap,v]  

            # update degree of neigbour nodes of v
            if v in neighbour_nodes:
                    # update degree
                    obj.deg -= 1                                
            
            deg = obj.deg
            # update density
            if n != 1:
                obj.den  = ((m - deg) / (n - 1) - (m / n))  
            else:
                obj.den = 0                                          

        #update pattern based sum contraint gains
        for idx in range(numberOfGraphs):
            snapObj = snap_dic[idx]
            pattern_dic = snapObj.pat_val

            for _pattern, val in pattern_dic.items():
                # print(_pattern)
                pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
        

        # update pattern based sum contraint gains
        # for idx in range(numberOfGraphs):
        #     snapObj = snap_dic[idx]
        #     pattern_dic = snapObj.pat_val
            
        #     if idx == max_snap:
        #         for _pattern, val in pattern_dic.items():
        #             # print(_pattern)
        #             pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )    
        #     else:
        #         i = max_snap
        #         _sum1 =(-lam*((I_old[i, idx] - 1)/U_old[i, idx] - (I_old[i, idx]/U_old[i, idx]))) + (lam*(((- 1)/(U[i, idx])) ))
        #         _sum2 = (-lam*(I_old[i, idx] /(U_old[i, idx] - 1)  - (I_old[i, idx]/U_old[i, idx]) )) + (lam*((I[i, idx]/((U[i, idx] - 1)*U[i, idx])) )) 
                
        #         for _pattern, val in pattern_dic.items():
        #             # print('idx {} - patt {}'.format(idx, _pattern))
        #             # print(' lst {} - {} : {}'.format(len(new_pat_lst), tuple(new_pat_lst[idx]),idx))
        #             if _pattern != tuple(new_pat_lst[idx]):
        #                 if _pattern[max_snap]:
        #                     pattern_dic[_pattern] += _sum1 
        #                     # print('val  : {} - {}'.format(pattern_dic[_pattern],lam*sum_con_pattern(I, U, idx, _pattern )) )
        #                 else:
        #                     pattern_dic[_pattern] += _sum2 
        #                     # print('val  : {} - {}'.format(pattern_dic[_pattern],lam*sum_con_pattern(I, U, idx, _pattern )) )
        
        

                        
    print("--- %s seconds ---" % (time.time() - start_time))
    print('count:')
    # print(cnt)
    # print(cnt_wrn)
    print(com_q(snapshots, set_dic, lam))
    return set_dic


def soft_2_org(snapshots, dcs, lam, nodes):
    list_of_sets = [set(copy.deepcopy(itm.nodes())) for itm in snapshots]
    numberOfGraphs = len(snapshots)
    set_dic = {}
    keys = range(numberOfGraphs)
    values = list_of_sets
    set_dic = dict(zip(keys, values))

    # set iterations
    list_of_sets_itr = [set(copy.deepcopy(itm.nodes())) for itm in snapshots]
    set_itr = {}
    keys = range(numberOfGraphs)
    values = list_of_sets_itr
    set_itr = dict(zip(keys, values))

    d = com_q(snapshots, set_dic, lam)
    print('current density: {}'.format(d))
    
    start_time = time.time()  
    
    # dictionary to store vertex objects
    v_count_dic = {}
    #dictionary to store snap objects
    snap_dic = {}
    
    # vertex counters initialization
    for idx, snap in enumerate(snapshots):

        I = []
        U = []
        
        for j in range(numberOfGraphs):
            if idx != j :
                I.append( len(set_dic[idx].intersection(set_dic[j])))
                # print('{} {}'.format(len(set_dic[idx].intersection(set_dic[j])), len(set_dic[idx].union(set_dic[j]))))
                U.append(len(set_dic[idx].union(set_dic[j])))
            else:
                I.append( len(set_dic[idx]))
                U.append( len(set_dic[idx]))
            
        snapObj =  Snapshot(idx)    
        m = len(snap.edges())
        n = len(snap.nodes())
        snapObj.m = m
        snapObj.n = n
        
        for i in snap:
            
            obj = Vertex((idx,i))
            
            # set vertex counters
            obj.deg = snap.degree(i)
            
            _pattern = []

            for ii, (k, lst1) in enumerate(set_dic.items()):
                if ii == idx:
                    _pattern.append(-1)
                elif i in lst1:
                    _pattern.append( 1)
                else:
                    _pattern.append( 0)
            obj.pattern = _pattern 
            # snapObj.AddPattern(tuple(_pattern), obj.key)
            snapObj.AddPattern(tuple(_pattern), i)
            
            v_count_dic[idx,i] = obj
        snap_dic[idx] = snapObj
        # snapObj.Display_Patterns()     
        # snapObj.Display_Patter_val()
    # print(v_count_dic[1,110])
         
    delta = 0
    maxx = 0

    # initial snapshots
    snapshots1 = copy.deepcopy(snapshots)
    # initialize I and U matrices
    I = {}
    U = {}
    for i in range(numberOfGraphs):
        for j in range(numberOfGraphs):
            if i < j :
                I[i, j] = len(set_dic[i].intersection(set_dic[j]))
                I[j, i] = I[i, j]  
                U[i, j] = len(set_dic[i].union(set_dic[j]))
                U[j, i] = U[i, j] 
                
    sum_cons_lst = []
    
    for idx1 in range(numberOfGraphs):
        _sum = 0
        for idx2 in range(numberOfGraphs):
            if idx1 != idx2:

                _sum += (lam*I[idx1,idx2]/U[idx1,idx2])
        sum_cons_lst.append(_sum  )  
        
    # initialization         
    for idx, snap_G in enumerate(snapshots1):

        if len(set_itr[idx]) > 0:

            S_i = set_itr[idx]
            
            snap = snap_dic[idx]
            
            m = snap.m
            n = snap.n

            for v in S_i:
                
                deg = v_count_dic[idx,v].deg
                
                if n != 1:
                    den = (m - deg) / (n - 1) - (m / n) 
                else:
                    den = 0

                v_count_dic[idx,v].den  = den
                            
    # for item in snap_dic.items():
    for idx in range(numberOfGraphs):
        snapObj = snap_dic[idx]
        pattern_dic = snapObj.pat_val
        
        for _pattern, val in pattern_dic.items():
            # print(_pattern)
            pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
            # print(pattern_dic[_pattern])
                   
    # iterate through all the nodes
    while not checkEmpty(set_itr):
        # find vertex with max total gain   (GREEDY DIFF (local))              
        _max = -math.inf
        max_snap = -1
        max_v = -1
        # print(numberOfGraphs)
        for _snap in range(numberOfGraphs):
                        
            snapObj = snap_dic[_snap] 
            patt_dic =  snapObj.pat_dic
            patt_val = snapObj.pat_val

            for pat in patt_val.keys():  
                    _max_pat = -math.inf
                    _max_pat_v= None
                    
                    #search for the vertex with max density component within a pattern
                    for v in patt_dic[pat]:
                        gain_v = v_count_dic[_snap,v].den
                        
                        if gain_v > _max_pat:
                            _max_pat = gain_v
                            _max_pat_v = v

                    tot_gain = v_count_dic[_snap,_max_pat_v].den  + patt_val[pat]
                        
                    if tot_gain > _max:
                        _max = tot_gain
                        
                        max_snap, max_v = _snap,_max_pat_v 
                        max_pattern = pat
                    
        # print('val: {} {}'.format(max_snap, max_v))        
            
        # current snapshots
        chg_snap = snapshots1[max_snap]
         # find neighbours of node v
        neighbour_nodes  = [n for n in chg_snap.neighbors(max_v)]
        
        diff =    _max
        # remove the vertex with the highest score
        #update SETS in each iteration
        set_itr[max_snap] = set_itr[max_snap] - {max_v}
        #update SNAPSHOT
        
        snapshots1[max_snap].remove_node(max_v)            
        
        delta += diff  
        # print('{} {}'.format(diff, obj.tot_gain))
        new_sets = copy.deepcopy(set_itr)
        
        # maintain GLOABL diff
        if delta  > maxx:
            set_dic = new_sets
            maxx = delta
       
        # U_old = copy.deepcopy(U)
        # I_old = copy.deepcopy(I)
        
        # # # update intersection and union matrices
        updateIU(I, U, set_dic, max_snap, max_v)
        
        # new_pat_lst = []
        
        #Update patterns of max_v from other snapshots
        for idx, snap_G in enumerate(snapshots1):
            
            if idx != max_snap:

                if max_v in snap_G:                    

                    snap_obj = snap_dic[idx]

                    pattern = tuple(v_count_dic[idx,max_v].pattern)
                    
                    # i = max_snap
                    # _sum = (-lam*(( - 1)/U_old[i, idx]  )) + (lam*((I[i, idx]/(U[i, idx]- 1)) - (I[i, idx]/U[i, idx])))
                    
                    new_pattern= snap_obj.updatePattern(max_v,idx,max_snap,pattern) 
                    v_count_dic[idx,max_v].pattern = new_pattern
                #     new_pat_lst.append(new_pattern)
                # else:
                #     new_pat_lst.append([0])
            # else:
            #      new_pat_lst.append([0])
                 
        #  max_snap snap Obj                           
        maxsnapObj = snap_dic[max_snap] 
        
        # update m and n counters after max_v is removed from max_snap
        maxsnapObj.m = maxsnapObj.m - v_count_dic[max_snap,max_v].deg
        maxsnapObj.n = maxsnapObj.n -  1
                
        #remove max_v from max_snap 
        max_pat_dic_lst =  maxsnapObj.pat_dic[max_pattern]
        max_pat_dic_lst.remove(max_v)
        
        if len(max_pat_dic_lst) == 0:
            maxsnapObj.pat_dic.pop(max_pattern)
            maxsnapObj.pat_val.pop(max_pattern)
        
        #     print('tot {} {}'.format(obj.tot_gain, diff))
        # remove max_v
        del v_count_dic[max_snap,max_v]
        
        # update new density gains in same snapshot        
        S_i = set_itr[max_snap]
                
        m = maxsnapObj.m
        n = maxsnapObj.n
        
        for v in S_i:
            obj  = v_count_dic[max_snap,v]  

            # update degree of neigbour nodes of v
            if v in neighbour_nodes:
                    # update degree
                    obj.deg -= 1                                
            
            deg = obj.deg
            # update density
            if n != 1:
                obj.den  = ((m - deg) / (n - 1) - (m / n))  
            else:
                obj.den = 0                                          

        #update pattern based sum contraint gains
        for idx in range(numberOfGraphs):
            snapObj = snap_dic[idx]
            pattern_dic = snapObj.pat_val

            for _pattern, val in pattern_dic.items():
                # print(_pattern)
                pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )
        

        # update pattern based sum contraint gains
        # for idx in range(numberOfGraphs):
        #     snapObj = snap_dic[idx]
        #     pattern_dic = snapObj.pat_val
            
        #     if idx == max_snap:
        #         for _pattern, val in pattern_dic.items():
        #             # print(_pattern)
        #             pattern_dic[_pattern] = lam*sum_con_pattern(I, U, idx, _pattern )    
        #     else:
        #         i = max_snap
        #         _sum1 =(-lam*((I_old[i, idx] - 1)/U_old[i, idx] - (I_old[i, idx]/U_old[i, idx]))) + (lam*(((- 1)/(U[i, idx])) ))
        #         _sum2 = (-lam*(I_old[i, idx] /(U_old[i, idx] - 1)  - (I_old[i, idx]/U_old[i, idx]) )) + (lam*((I[i, idx]/((U[i, idx] - 1)*U[i, idx])) )) 
                
        #         for _pattern, val in pattern_dic.items():
        #             # print('idx {} - patt {}'.format(idx, _pattern))
        #             # print(' lst {} - {} : {}'.format(len(new_pat_lst), tuple(new_pat_lst[idx]),idx))
        #             if _pattern != tuple(new_pat_lst[idx]):
        #                 if _pattern[max_snap]:
        #                     pattern_dic[_pattern] += _sum1 
        #                     # print('val  : {} - {}'.format(pattern_dic[_pattern],lam*sum_con_pattern(I, U, idx, _pattern )) )
        #                 else:
        #                     pattern_dic[_pattern] += _sum2 
        #                     # print('val  : {} - {}'.format(pattern_dic[_pattern],lam*sum_con_pattern(I, U, idx, _pattern )) )
        
        

                        
    print("--- %s seconds ---" % (time.time() - start_time))
    print('count:')
    # print(cnt)
    # print(cnt_wrn)
    print(com_q(snapshots, set_dic, lam))
    return set_dic