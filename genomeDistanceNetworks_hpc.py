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
import gzip


seed = np.random.RandomState(seed=3)
distance_dict = defaultdict(dict)


## load distance dictionary and labels

#
# dataFile = 'distances_filtered.csv'
# projectionFile = 'proj.pkl'
# distMatrixFile = 'distMatrix_norm.pkl'

os.chdir('/home/u1474301/fastani/')
ref_genomes = load(open('ref_acc.pkl','rb'))
ref_genomes = set(ref_genomes)
#pairs_file = '/Users/emzodls/Documents/pairs_9868_cutoff.pkl'
#species_file = '/Users/emzodls/Dropbox/Lab/Warwick/distance-tests/species.pkl'

distFile = 'NCBI_Prok-matrix.txt.gz'
species = set()
pairs = []
ref_pairs = []
max_similarity = 0
cutoff = 90
logfile  = 'file_log.out'
for line in gzip.open(distFile,'rt'):
    lineParse = line.split()
    species1 = '_'.join(lineParse[0].split('.')[0].split('_')[-2:])
    species2 = '_'.join(lineParse[1].split('.')[0].split('_')[-2:])
    species.add(species1)
    species.add(species2)
    similarity = float(lineParse[-1])
    if species1 and species2 in ref_genomes and species1 != species2:
        ref_pairs.append((species1,species2,similarity))
        if similarity >= max_similarity:
            max_similarity = similarity
    if similarity >= cutoff:
        pairs.append((species1,species2,similarity))

dump(ref_pairs,open('ref_pairs.pkl','wb'))
print('Finished Reading File')
print(len(ref_pairs))
print(max_similarity)
species = sorted(list(species))
species2Idx = dict(zip(species,range(len(species))))
idx2species = {v:k for k,v in species2Idx.items()}
with open(logfile,'w') as outfile:
    outfile.write('Finished Reading File\n')
    outfile.write('Found {} reference genome comparisons, with a maximum similarity of {}\n'.format(len(ref_pairs),max_similarity))
simDict = dict()
for (species1,species2,similarity) in pairs:
    if species1 != species2 and (similarity >= max_similarity):
        species_order = tuple(sorted([species1,species2]))
        stored_similarity  = simDict.setdefault(species_order,0)
        if similarity >= stored_similarity:
            simDict[species_order] = similarity

dump(simDict,open('simDict.pkl','wb'))
## generate network

g = igraph.Graph()
print('Found {} different assemblies, converting to nodes'.format(len(species)))
for genome in species:
    g.add_vertex(genome)
totalEdges = len(simDict)
# for edgenum,(species1,species2) in enumerate(simDict.keys()):
#     print('Added {} edge of {}'.format(edgenum,totalEdges))
#     g.add_edge(species1,species2, weight=simDict[(species1,species2)])
edgesAndWeights = [(k,v) for k,v in simDict.items()]
g.add_edges([x for x,y in edgesAndWeights])
g.es["weight"] = [y for x,y in edgesAndWeights]
## get maximal cliques
#g.maximal_cliques(file=open('/Users/emzodls/Documents/cliques.txt','w'))
print('Added {} edges'.format(len(simDict)))
with open(logfile,'a') as outfile:
    outfile.write('Added {} edges\n'.format(len(simDict)))
### Group with label propagation
subgraphs = g.community_multilevel()
membership = subgraphs.membership
with open('node_membership_multilevel.csv','w') as outfile:
    for vertex,membership in zip(g.vs,membership):
        outfile.write('{},{}\n'.format(vertex["name"],membership))

with open(logfile,'a') as outfile:
    outfile.write('Grouped into {} groups\n'.format(len(set(subgraphs.membership))))
# ### Group with walkthrough
# dendrogram = g.community_walktrap()
# clusters_walktrap = dendrogram.as_clustering()
# membership_walktrap = clusters_walktrap.membership
# print('Finding neighborhoods with walktrap')
# with open('/Users/emzodls/Documents/node_membership_walktrap.csv','w') as outfile:
#     for vertex,membership in zip(g.vs,membership_walktrap):
#         outfile.write('{},{}\n'.format(vertex["name"],membership))
