import os
import sys
sys.path.append(os.path.abspath(os.path.join("collect_tf_motifs.json" ,"../..")))
import json
from Bio import motifs
from soa_operon_cluster import OperonCluster
from soa_motif_analysis import find_pattern, calculate_motif_distance,calc_kld_distance, calc_euclidean, write_to_csv

input_file = "test_operons.txt"
dist_calc = calc_euclidean
#dist_calc = calc_kld_distance

# Create operon clusters consisting of only one motif each (every motif consists of only one site)
# as specified by the input file
def get_clusters():
    f = open(input_file, "r")
    lines = f.readlines()

    operon_clusters = []
    for line in lines:
        line = line.split(':')
        array = []
        array.append(line[1].replace('\n',''))
        motif = motifs.create(array)
        motifs_with_repeats = []
        if (find_pattern(motif)):
            motifs_with_repeats.append(motif)

        new_cluster = OperonCluster(cluster_id = line[0])
        new_cluster.motifs = motifs_with_repeats
        operon_clusters.append(new_cluster)

    return operon_clusters

if __name__ == "__main__":
    clusters = get_clusters()
    write_to_csv(input_file.split('.')[0], clusters, dist_calc)
