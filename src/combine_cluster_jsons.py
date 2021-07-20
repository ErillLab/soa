'''
@author ichaudr
7/12/21
Objective: Combine the outputs of different runs.
'''

import json
import os

#Directories
to_read_in = ['../output/complete_clusters_85pcent/', '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/archive/complete_clusters_06162021/complete_clusters/']

#Final list of clusters that will be outputted
final_out = {}

for path in to_read_in:
    files_in_path = list(os.listdir(path))
    for file in files_in_path:

        if file.split('.')[-1] != 'json':
            continue

        cluster_info = dict(json.load(open(path+file, 'r')))

        if cluster_info['cluster_id'] in final_out.keys():
            final_out[cluster_info['cluster_id']]['motifs'].extend(cluster_info['motifs'])
        else:
            final_out[cluster_info['cluster_id']] = cluster_info

#Output directory
output_dir = '../output/combined_85_75_07122021/'
for clster in final_out.keys():
    with open(output_dir + clster + '.json', 'w+') as file:
        json.dump(final_out[clster], file)