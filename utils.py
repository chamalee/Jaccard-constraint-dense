
from datetime import datetime, timedelta
import time
from gurobipy import *
import networkx as nx


def readFile(filepath):

    nodes = set()
    edges = set()
    time_edges = [] 

    with open(filepath,'r') as fd:
        for l in fd.readlines():
            
            l = l.strip()            
            items = l.split(' ')            
            tstamp = ' '.join(items[0:2])
            tstamp = tstamp[1:-1]

            tstamp = int(time.mktime(datetime.strptime(tstamp, '%Y-%m-%d %H:%M:%S').date().timetuple()))

            t = items[2:4]
            t = list(map(int,t))

            if t[0] == t[1]:            
                continue
            t.sort(); 

            time_edges.append(( t[0], t[1],tstamp ))

            nodes.add(t[0])
            nodes.add(t[1])
            edges.add(tuple([t[0],t[1]]))

    fd.close()
    return time_edges, nodes, edges
    
