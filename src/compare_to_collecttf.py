'''
@author ichaudr

Comparing the motifs pulled by the pipeline to the motifs that are determined on collecttf.

'''

from Bio.motifs import Motif, Instances
import csv
import json
import os


## Helper functions
def get_cluster_motifs(cluster_id):
    '''
    Parse out the motif for a specific cluster.

    Parameters
    ----------
    cluster_id: str
        The cluster of interest.
    path_to_saved_clusters: str
        The path to the directory that is holding the saved JSONs for the saved clusters
    
    Returns
    -------
    motifs: str
        A string that is made of the concatenated consensus sequences of the motifs a part of the cluster. 
    '''
    path_to_saved_clusters = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/complete_clusters/'
    motifs = []

    for file in os.listdir(path_to_saved_clusters):
        if file.split('.')[0] == cluster_id:
            file_reader = json.load(open(path_to_saved_clusters+file, 'r'))

            motifs_from_file = file_reader['motifs']

            for m in motifs_from_file:
                temp_motif = Motif(instances=Instances(m))
                motifs.append(str(temp_motif.consensus))
    return motifs




def get_cluster(locus_tag):
    '''
    Find the cluster a specific gene belongs to.
    '''

    #Path to the cluster info json
    cluster_info_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/src/edge_file_annotator/cluster_defs.json'

    with open(cluster_info_path, 'r') as file:
        clusters = json.load(file)

        for c in clusters:
            for g in clusters[c].values():
                if g['locus_tag'] == locus_tag or g['old_locus_tag'] == locus_tag:
                    return c
    
    return None





#Path to collecttf data
#collecttf_data_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/EdgeFileAnalysis/collectf-export-tsv-raw.csv'
collecttf_data_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/misc/collecttf/vibrio_tf_info.csv'

#Parse out tf, motifs, and regulated genes
tf_data = []

with open(collecttf_data_path, 'r') as file:
    reader = list(csv.reader(file))
    for row in reader[1:]:
        tf = row[0]
        #tf_acc = row[1]
        #genome_acc = row[2]
        tf_motif_seq = row[1]
        tf_genes = list(row[2].replace("'",'').replace('[','').replace(']','').split(','))
        tf_data.append({
            'tf':tf,
            #'tf_acc':tf_acc,
            #'genome_acc':genome_acc,
            'tf_motif_seq':tf_motif_seq,
            'tf_genes':tf_genes,
            'predicted_motif_seqs':[]
        })

for tf in tf_data:
    genes = tf['tf_genes']
    predicted_motifs = []
    for g in genes:
        if get_cluster(g) == None:
            continue
        predicted_motifs.extend(get_cluster_motifs(get_cluster(g)))
    tf['predicted_motif_seqs'] = ', '.join(predicted_motifs)


#Write to output
ouput_file = 'motif_compare2.csv'

with open(ouput_file, 'a') as file:
    writer = csv.writer(file)
    writer.writerow(list(tf_data[0].keys()))

    for tf in tf_data:
        if len(tf['predicted_motif_seqs']) == 0:
            continue
        writer.writerow(list(tf.values()))
    