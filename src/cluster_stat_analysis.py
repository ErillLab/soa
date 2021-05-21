'''
@author ichaudr

This script will report a set of statistics for all of the clusters that were saved. Specifically, it will report:

- Average number of motifs per cluster
- Max number of motifs per cluster
- Min number of motifs per cluster
- Median number of motifs per cluster

- Average IC of all the motifs
- Max IC of all of the motifs
- Min IC of all of the motifs
- Median IC of all of the motifs

'''

from Bio.motifs import Motif, Instances
import json
import csv
import os
import statistics

#Path to the cluster JSONs
path_to_clusters = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/complete_clusters/'

#A dictionary holding the number of motifs per cluster
num_motifs_per_cluster = {}

#A dictionary holding all of the motifs and their IC per cluster
ic_per_motif_per_cluster = {}

#Get the counts to fill the above dictionaries
for file in os.listdir(path_to_clusters):

    if file[-4:] != 'json':
        continue
    
    cluster_id = file[:-5]
    cluster_reader = json.load(open(path_to_clusters+file, 'r'))

    #Add info about the number of motifs
    num_motifs_per_cluster[cluster_id] = len(cluster_reader['motifs'])

    #Add info about the IC per motif
    ic_per_motif_per_cluster[cluster_id] = {}
    for m in range(len(cluster_reader['motifs'])):
        temp_motif = Motif(instances=Instances(cluster_reader['motifs'][m]))
        temp_motif_seq = str(temp_motif.consensus)
        temp_motif_ic = temp_motif.pssm.mean()

        ic_per_motif_per_cluster[cluster_id][temp_motif_seq] = temp_motif_ic


#All stats to be calculated:
avg_num_motifs_per_cluster = statistics.mean(list(num_motifs_per_cluster.values()))
stdv_num_motifs_per_cluster = statistics.stdev(list(num_motifs_per_cluster.values()))
max_num_motifs_per_cluster = max(list(num_motifs_per_cluster.values()))
min_num_motifs_per_cluster = min(list(num_motifs_per_cluster.values()))
median_num_motifs_per_cluster = statistics.median(list(num_motifs_per_cluster.values()))
mode_num_motifs_per_cluster = statistics.mode(list(num_motifs_per_cluster.values()))

all_ics = []
for c in ic_per_motif_per_cluster.keys():
    for m in ic_per_motif_per_cluster[c].keys():
        all_ics.append(ic_per_motif_per_cluster[c][m])

avg_ic = statistics.mean(all_ics)
stdv_ic = statistics.stdev(all_ics)
max_ic = max(all_ics)
min_ic = min(all_ics)
median_ic = statistics.median(all_ics)
mode_ic = statistics.mode(all_ics)

print('-'*10)
print('Statistics for the number of motifs per cluster (n=', len(list(num_motifs_per_cluster.values())),'):')
print('\tAverage number of motifs per cluster: ', avg_num_motifs_per_cluster, ' ± ', stdv_num_motifs_per_cluster)
print('\tMax number of motifs per cluster: ', max_num_motifs_per_cluster)
print('\tMin number of motifs per cluster: ', min_num_motifs_per_cluster)
print('\tMedian number of motifs per cluster: ', median_num_motifs_per_cluster)
print('\tMode of number of motifs per cluster: ', mode_num_motifs_per_cluster)
print('-' * 10)
print('Statistics for the IC content of motifs (n=', len(all_ics), '):')
print('\tAverage IC: ', avg_ic, ' ± ', stdv_ic)
print('\tMax IC: ', max_ic)
print('\tMin IC: ', min_ic)
print('\tMedian IC: ', median_ic)
print('\tMode IC: ', mode_ic)
print('-'*10)

with open('./motif_counts.csv', 'w') as file:
    filew = csv.writer(file)
    for c in list(num_motifs_per_cluster.values()):
        filew.writerow([c])

with open('./motif_ic.csv', 'w') as file:
    filew = csv.writer(file)
    for ic in all_ics:
        filew.writerow([ic])
