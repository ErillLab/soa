'''
An accessory script to calculate the distances from the saved clusters.

@autho ichaudr
04.20.2021
'''

from soa_motif_analysis import *
from soa_operon_cluster import *
import os

cluster_dir_path = '/home/issac/Desktop/complete_clusters_04202021/'

#Holds all of the operon clusters
all_clusters =[]

for file in os.listdir(cluster_dir_path):
	if file.split('.')[-1] != 'json':
		continue
	print('Adding: ', str(file))
	c = OperonCluster()
	c.load_from_json(cluster_dir_path + file)
	all_clusters.append(c)

write_to_csv(filename='complete_clusters_04202021_test', operon_clusters=all_clusters, distance_calc=calc_euclidean)
#write_to_csv(filename='complete_clusters_04202021_kld_only', operon_clusters=all_clusters, distance_calc=calc_kld_distance)