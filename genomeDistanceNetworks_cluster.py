# -*- coding: utf-8 -*-
"""
Created on Thu 11 Feb 2016

@author: emzodls
"""


from collections import defaultdict
import numpy as np
import pandas as pd
import sys,os
from pickle import dump,load
import igraph

seed = np.random.RandomState(seed=3)
distance_dict = defaultdict(dict)


## load distance dictionary and labels

#
# dataFile = 'distances_filtered.csv'
# projectionFile = 'proj.pkl'
# distMatrixFile = 'distMatrix_norm.pkl'

os.chdir('/home/u1474301/ncbi_prok_derep/')

pairs_file = 'pairs_9868_cutoff.pkl'
species_file = 'species.pkl'


# pairs = load(open(pairs_file,'rb'))
elements = load(open(species_file,'rb'))

elements = sorted(list(elements))
elements2Idx = dict(zip(elements,range(len(elements))))
idx2elements = {v:k for k,v in elements2Idx.items()}
simDict = load(open('simDict_9868_cutoff.pkl','rb'))

## generate network
print('Found {} different assemblies, converting to nodes'.format(len(elements)))
g = igraph.Graph()

for species in elements:
    g.add_vertex(species)

for (species1,species2) in simDict.keys():
    g.add_edge(species1,species2, weight=simDict[(species1,species2)])
print('Added {} edges, converting to nodes'.format(len(simDict)))
## get maximal cliques
g.maximal_cliques(file=open('/Users/emzodls/Documents/cliques.txt'))

### Group with walkthrough
print('Finding neighborhoods with walktrap')
dendrogram = g.community_walktrap(weights="weight")
clusters_walktrap = dendrogram.as_clustering()
membership_walktrap = clusters_walktrap.membership
print('Writing Neighborhood file')
with open('node_membership.txt','w') as outfile:
    for vertex,membership in zip(g.vs,membership_walktrap):
        outfile.write('{},{}\n'.format(vertex["name"],membership))


