'''
@author ichaudr

A testing module to test other modules.
'''
import os
import sys
sys.path.append(os.path.abspath(os.path.join("testing_motif_distances.py" ,"../../")))
import soa_sim_filter
from soa_motif_analysis import  *
from Bio.motifs import Motif, Instances


test_meme_output_1 = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/meme_bin/CLSTR1708_meme_out/'
test_meme_output_2 = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/meme_bin/CLSTR771_meme_out/'

records_1 = get_motifs(test_meme_output_1, 100)
records_2 = get_motifs(test_meme_output_2, 100)

#print(len(records_1))

#print(len(records_2))

#motif_a = records_1[2]
#motif_b = records_2[0]

motif_a = records_1[0]
motif_b = records_2[0]

#motif_a = Motif(instances=Instances(['AAAAAAAAAAAAAAAAAAAA']))
#motif_b = Motif(instances=Instances(['AAAAAAAAAAAAAAAAAAAG']))
#motif_b = Motif(instances=Instances(['TGCTCGATCGATCGATTGCTA']))

#print(motif_a.pwm)
#print(motif_b.pwm)

print(calculate_motif_distance(motif_a, motif_b, padded=True, distance_function=calc_kld_distance, add_psuedocounts=False, scaling_factor=10))

val = .25
#val = (len(str(motif_b.instances).split('\n')[0]) + len(str(motif_a.instances).split('\n')[0]) / 2)/(len(motif_b.instances) + len(motif_a.instances) / 2)
test_dict = dict(A=val, C=val, G=val, T=val)

#motif_a.pseudocounts = test_dict
#motif_b.pseudocounts = test_dict

print(calculate_motif_distance(motif_a, motif_b, padded=True, distance_function=calc_kld_distance, add_psuedocounts=True, scaling_factor=10))


#print(motif_a.pwm)
#print(motif_b.pwm)
