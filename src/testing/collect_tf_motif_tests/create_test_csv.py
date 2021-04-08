'''
Generates Gephi csv files from test operons consisting of motifs from CollectTF as specified by an input file using 
the specified motif distance calculation method (either euclidean or kld)

Commands to run:
    python create_test_csv.py sample_operons+motifs.txt e/k
    * e = euclidean distance calculation between motifs
    * k = kld distance calculation between motifs
    * Defaults to hardcoded values if no argumentss are passed in

Input file format (see sample_operons+motifs.txt as a reference):
    OPRN [insert name]
    motif one (tf_instance name from Collect TF)
    motif two (tf_instance name from Collect TF)
    OPRN [insert name]
    motif one (tf_instance name from Collect TF)
    motif two (tf_instance name from Collect TF)

CSV files with the following operon-distance calculations are created:
    a. average of the distances between ALL motifs in both operons
    b. minimum of the distances between ALL motifs in both operons
    c. maximum of the distances between ALL motifs in both operons
    d. median of the distances between ALL motifs in both operons
    e. noise reduction distance calculation
        Let the operon with fewer motifs be Operon A, and the other be Operon B. For each motif in Operon A, 
        find the motif in Operon B for which the motif in Operon A has the smallest distance. Then, take
        the average of these minimum distances to calculate the distance between the operons
'''

import os
import sys
sys.path.append(os.path.abspath(os.path.join("collect_tf_motifs.json" ,"../../..")))
import json
from Bio import motifs
from soa_operon_cluster import OperonCluster
from soa_motif_analysis import find_pattern, calculate_motif_distance,calc_kld_distance, calc_euclidean, write_to_csv

# Set default file (in correct format as described above) containing operons and their motifs from Collect TF
operon_file = 'sample_operons+motifs.txt'
# Set defaul motif distance calculation - either 'calc_euclidean' or 'calc_kld_distance'
#distance_calc = calc_euclidean
distance_calc = calc_kld_distance

def get_clusters_from_file():
    '''
    Creates operon cluster objects from specified file containing operons and their motifs from Collect TF

    Parameters
    ----------
    None

    Returns
    -------
    operon_clusters:[operon_cluster] 
        An array holding operon_cluster objects for each operon in the specifed file
    '''

    # open json file with all Collect TF motifs
    with open('collect_tf_motifs.json') as f:
        collect_tf_motifs = json.load(f)

    # create operon cluster objects and append to operon_clusters
    operon_clusters = []
    with open(operon_file) as fp:
        line = fp.readline().strip()
        new_cluster = OperonCluster(cluster_id = line)

        while line:
            if (line[0:4].strip() == 'OPRN'):
                operon_clusters.append(new_cluster)
                new_cluster = OperonCluster(cluster_id = line)
            else:
                # Search for instance in collect_tf_motifs
                motif = create_motifs(line, collect_tf_motifs)
                # If tf_instance was found, add to operon_cluster's motifs
                if (motif != -1):
                    new_cluster.motifs.append(motif)
            line = fp.readline().strip()

        operon_clusters.append(new_cluster)
        operon_clusters.pop(0)
    return operon_clusters


def create_motifs(tf_instance_id, collect_tf_motifs):
    '''
    Searches json for the motif instance passed into the function and constructs motif object
    from the instance's aligned binding site. 

    Parameters
    ----------
    tf_instance_id: String

    Returns
    -------
    motif: motifs object
        Constructed motif object from aligned binding sites, name is set to tf instance id
    -1: Int
        Returns -1 when the tf_instance was not found
    '''
    entry = list(filter(lambda motif: motif['tf_instance'] == tf_instance_id, collect_tf_motifs['all_motifs']))
    if (len(entry) > 0):
        motif = motifs.create(entry[0]['aligned_binding_sites'])
        motif.name = entry[0]['tf_instance']
        return motif
    else:
        return -1

if __name__ == "__main__":
    # Get command line arguments if passed in, otherwise use hardcoded values
    if (len(sys.argv) > 1):
        # Get and set input file name
        operon_file = sys.argv[1]
        # Get and set distance_calc function
        if (sys.argv[2][0] == 'e'):
            print(sys.argv[2][0])
            distance_calc = calc_euclidean
        elif (sys.argv[2][0] == 'k'):
            print(sys.argv[2][0])
            distance_calc = calc_kld_distance

    # Create operon_cluster objects from input file
    clusters = get_clusters_from_file()

    # Print operon name ... motifs
    for cluster in clusters:
        for motif in cluster.motifs:
            print(cluster.cluster_id,'...', motif.name)

    output_file = operon_file.split('.')[0]
    write_to_csv(output_file, clusters, distance_calc)