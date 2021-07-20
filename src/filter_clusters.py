'''
@author ichaudr
July 17 2021

Objective: Takes in clusters with unfiltereing promoters and filters them at a given threshold similarity level. The motifs are discovered via MEME after the filtering and the clusters are re-written as jsons.

'''

import json
import os

from soa_sim_filter import sim_filter
from soa_motif_analysis import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tqdm import tqdm

#threshold value
THRESHOLD = 0.75

#path to input clusters, these are unfiltered
unfiltered_clusters_path = 'complete_clusters_no_filtering/'

#read in unfiltered clusters
unfiltered_clusters = []
for file in os.listdir(unfiltered_clusters_path):
    if file.split('.')[-1] != 'json':
        continue
    unfiltered_clusters.append(dict(json.load(open(unfiltered_clusters_path + file, 'r'))))

#unfiltered_clusters = [dict(json.load(open('complete_clusters_no_filtering/CLSTR643.json', 'r')))]

#iterate through each of the clusters
for cluster in tqdm(unfiltered_clusters):

    tqdm.write(cluster['cluster_id'])

    #parse out all the promoters for each of the operons in the cluster
    cluster_promoters = [operon['promoter'] for operon in cluster['operons']]

    #filter the promoters at the set threshold value
    filtered_promoters_seqs =  sim_filter(sequences=cluster_promoters, threshold_percent_id=THRESHOLD)
    filtered_promoters = []

    i = 0
    for prom in filtered_promoters_seqs:
        filtered_promoters.append(SeqRecord(Seq(prom), id=cluster['cluster_id']+'_'+'prom'+str(i)))
        i+=1

    cluster['filtered_promoters'] = [str(f.seq) for f in filtered_promoters]

    #write promoters to temp fasta file
    SeqIO.write(filtered_promoters, 'temp.fasta', 'fasta')

    if len(filtered_promoters) < 5:
        os.system('rm temp.fasta')
        continue

    #meme discovery on the filtered promoters
    run_meme('temp.fasta', 'meme_output_for_cluster_filtering', width_min=10, width_max=20)
    os.system('rm temp.fasta')
    cluster_motifs = get_motifs('meme_output_for_cluster_filtering/', e_val_threshold=10E-3)

    #add the discovered motifs to the cluster dictionary
    cluster['motifs'].extend([str(m.instances).split('\n')[:-1] for m in cluster_motifs if (len(str(m.instances).split('\n')) > 5 and find_pattern(target_motif=m))])

    if len(cluster['motifs']) == 0:
       continue

    #write cluster to output directory
    with open('filtered_clusters_output/' + cluster['cluster_id'] + '.json', 'w+') as file:
        print('\tWriting ', cluster['cluster_id'])
        json.dump(cluster, file)
    
