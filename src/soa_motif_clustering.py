'''
@author ichaudr
Objective is to cluster motifs that are most similar together into bins. This information will be used to annotated the edges in the final SOA network. 

Approach:
    1. Import all of the motifs from the saved clusters. 
    2. Bin the motifs.
    3. Go through each of the clusters and annotate in a seperate file the bins the clusters' motifs are a part of.
'''

from Bio.motifs import Motif, Instances
from soa_motif_analysis import *

from tqdm import tqdm

import os
import json


### STEP 1: Import all of the motifs from the saved clusters

#Path to the operon_clusters
operon_clusters_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/final_soa_and_analysis/pipeline_output/complete_clusters/'

#Holds all the motifs from all the clusters
all_motifs = []

#Holds motifs per cluster
motifs_per_cluster = {}

#Import the motifs
for file in tqdm(list(os.listdir(operon_clusters_path)), desc='Reading in motifs from saved clusters'):
    if file.split('.')[-1] != 'json':
        continue
    
    cluster_info = dict(json.load(open(operon_clusters_path + file, 'r')))
    motifs_per_cluster[cluster_info['cluster_id']] = []

    for motif in cluster_info['motifs']:
        motif_obj = Motif(instances=Instances(motif))
        all_motifs.append(motif_obj)
        motifs_per_cluster[cluster_info['cluster_id']].append(motif_obj)

print('-'*5, ' Total motifs: ', len(all_motifs))

### STEP 2: Bin the motifs based off some threshold weight. 

#Holds all motif bins
print('Binning motifs...')
sim_threshold = 0.22
motif_bins = {}

bin_id = 0
while len(all_motifs) > 0:
    print('-'*5, ' Making bin ', bin_id, ' ... ', len(all_motifs), ' motifs remaining...')
    curr_motif = all_motifs.pop(0)
    motif_bins[bin_id] = [curr_motif]

    temp_motifs = []

    for other_motif in all_motifs:
        if calculate_motif_distance(curr_motif, other_motif, distance_function=calc_kld_distance) >= sim_threshold:
            motif_bins[bin_id].append(other_motif)
        else:
            temp_motifs.append(other_motif)
    
    all_motifs = temp_motifs
    temp_motifs = []
    bin_id += 1

print('Made ', len(motif_bins), ' bins.')


### STEP 3: Write output with the bin IDs present in each of the initial operon clusters.

#Holds the bin ids for each cluster...will be written to ouput file
motif_ids_per_cluster = {}

for cluster in motifs_per_cluster.keys():
    motif_ids_per_cluster[cluster] = []

    for motif in motifs_per_cluster[cluster]:
        for bin_id in motif_bins.keys():
            if motif in motif_bins[bin_id]:
                motif_ids_per_cluster[cluster].append(str(bin_id))

print('Output: ')
output = open('motif_clustering_out_' + str(sim_threshold).replace('.', '_') + '.txt', 'w+')
for cluster in motif_ids_per_cluster.keys():
    to_write = cluster + ':' + ','.join(motif_ids_per_cluster[cluster])
    output.write(to_write + '\n')
    print(to_write)
    




