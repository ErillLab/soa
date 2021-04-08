'''
@author ichaudr

Testing the Cluster.load_from_json function.
'''

import os
import sys
sys.path.append(os.path.abspath(os.path.join("test_cluster_loading.py" ,"../../")))
import soa_sim_filter
from soa_motif_analysis import  *
from Bio.motifs import Motif, Instances

c = OperonCluster()
c.load_from_json(file_path='/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/complete_clusters_04072021/CLSTR1602.json')
print(c)
