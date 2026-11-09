"""
Adaptation of SAND for transcripts: identification of local structures in spatial data

Inputs  1) X: Matrix where columns are x,y spatial coordinates and rows are transcripts
        2) r: distance cutoff for the construction of a radius-based graph
        
Returns 1) Spatial clusters of transcripts in X
"""

import pandas as pd
import pickle
import anndata as ad
from scipy.cluster.hierarchy import DisjointSet
from sklearn.neighbors import NearestNeighbors
import scipy

def runSAND_transcript(X, r):
    filtered_neighs = makeneighborgraph(X, r)
    return joinsets(filtered_neighs)
    
def makeneighborgraph(X, r): # X = x,y coords; r = radius
    tree = scipy.spatial.KDTree(X)
    k = tree.query_ball_tree(tree, r) 

    filtered_neighs = dict()
    for i,neighs in enumerate(k):
        filtered_neighs[i] = set(neighs) 
        
    return filtered_neighs

def joinsets(filtered_neighs): # call union find
    edges = set()
    nodes = set()
    for k in filtered_neighs:
        nodes.add(k)
        for v in filtered_neighs[k]:
            edges.add((k,v))
            
    disjoint_set = DisjointSet(nodes)
    for a,b in edges:
        disjoint_set.merge(a, b)

    return disjoint_set.subsets()