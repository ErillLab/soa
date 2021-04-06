from Bio import motifs
from soa_operon_cluster import OperonCluster
from soa_motif_analysis2 import find_pattern, calculate_motif_distance,calc_kld_distance
import os
import statistics
operon_clusters = []

filename = "test_motifs_2.txt"

def write_to_edge_file(file_handle, cluster_one, cluster_two, value):
    file_handle.write(cluster_one)
    file_handle.write(',')
    file_handle.write(cluster_two)
    file_handle.write(',')
    file_handle.write(value)
    file_handle.write('\n')

f = open(filename, "r")
lines = f.readlines()

for line in lines:
    line = line.split(':')
    array = []
    array.append(line[1].replace('\n',''))
    motif = motifs.create(array)
    motifs_with_repeats = []
    if (find_pattern(motif.instances)):
        motifs_with_repeats.append(motif)
        #print(motif.instances)

    new_cluster = OperonCluster(cluster_id = line[0])
    new_cluster.motifs = motifs_with_repeats
    operon_clusters.append(new_cluster)

for cluster in operon_clusters:
    for motif in cluster.motifs:
        print(motif.instances)
# Create file objects for various distance calculation methods (for edge csv files)
'''
edges_output_file_average = open("edges_average.csv", "w")
edges_output_file_average.write("Source,Target,Weight\n")

edges_output_file_median = open("edges_median.csv", "w")
edges_output_file_median.write("Source,Target,Weight\n")

edges_output_file_minimum = open("edges_minimum.csv", "w")
edges_output_file_minimum.write("Source,Target,Weight\n")

edges_output_file_maximum = open("edges_maximum.csv", "w")
edges_output_file_maximum.write("Source,Target,Weight\n")
'''
edges_output_file_noise_reduction = open(filename +"_edges_noise_reduction.csv", "w")
edges_output_file_noise_reduction.write("Source,Target,Weight\n")

# Create file oject for node csv file
nodes_output_file = open(filename + "_nodes.csv", "w")
nodes_output_file.write("Id,Label\n")

for i in range(0,len(operon_clusters)):
    # Create a node for every cluster in nodes.csv
    nodes_output_file.write(operon_clusters[i].cluster_id)
    nodes_output_file.write(',')
    nodes_output_file.write(operon_clusters[i].cluster_id)
    nodes_output_file.write("\n")

    # Compare the motifs of each cluster with the motifs of the other cluster
    for j in range(i,len(operon_clusters)):
        # Don't compare a cluster with itself
        if (i != j):
            # Hold all pairwise distances and minimum distances
            pairwise_distances = []
            minimum_distances = []

            # Identify the operon with fewer motifs
            # Always search for the minimum distance between a motif from the operon with the 
            # fewer amount of motifs and all of the motifs of the operon with the larger amount of motifs (noise reduction)
            if(len(operon_clusters[i].motifs) <= len(operon_clusters[j].motifs)):
                operon_fewer_motifs = operon_clusters[i].motifs
                operon_more_motifs = operon_clusters[j].motifs
            else:
                operon_fewer_motifs = operon_clusters[j].motifs
                operon_more_motifs = operon_clusters[i].motifs

            # Compare all of the motifs in cluster one with the motifs of cluster two
            for k in range(0, len(operon_fewer_motifs)):
                # Set initial minimum value if a comparison can be made
                if (len(operon_more_motifs) > 0):
                    min_val = calculate_motif_distance(operon_fewer_motifs[k],operon_more_motifs[0])
                
                for l in range(0, len(operon_more_motifs)):
                    # calculate kld/euclidean distance
                    distance_between_motifs = calculate_motif_distance(operon_fewer_motifs[k],operon_more_motifs[l], distance_function=calc_kld_distance)
                    print(distance_between_motifs)
                    # store all calculated distances
                    pairwise_distances.append(distance_between_motifs)

                    # If found new minimum, set min_val as new minimum value
                    if (distance_between_motifs < min_val):
                        min_val = distance_between_motifs

                    # If last comparison, add final min_val to minimum_distances array
                    if (l == len(operon_more_motifs)-1):
                        minimum_distances.append(min_val)
                    
            # If minimum distances for noise reduction were calculated
            if (len(minimum_distances) > 0):
                write_to_edge_file(edges_output_file_noise_reduction,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(1-sum(minimum_distances)/len(minimum_distances)))
                #print(i, "/" , len(operon_clusters), "--", operon_clusters[i].cluster_id, "|",operon_clusters[j].cluster_id,str(sum(minimum_distances)/len(minimum_distances)))
            # if pairwise_distances were calculated, write average, median, min, and max values to csv filie
            '''
            if len(pairwise_distances) > 0:
                write_to_edge_file(edges_output_file_average,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(sum(pairwise_distances)/len(pairwise_distances)))
                write_to_edge_file(edges_output_file_median,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(statistics.median(pairwise_distances)))
                write_to_edge_file(edges_output_file_minimum,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(min(pairwise_distances)))
                write_to_edge_file(edges_output_file_maximum,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(max(pairwise_distances)))
                print(i, "/" , len(operon_clusters), "--", operon_clusters[i].cluster_id, "|",operon_clusters[j].cluster_id)
            '''
# close all output files
'''
edges_output_file_average.close()
edges_output_file_median.close()
edges_output_file_maximum.close()
edges_output_file_minimum.close()
'''
edges_output_file_noise_reduction.close()
nodes_output_file.close()

edges_output_file_noise_reduction.close()
nodes_output_file.close()