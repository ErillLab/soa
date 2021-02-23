'''
@author ichaudr

A testing module to test other modules.
'''

import soa_sim_filter
from soa_motif_analysis import get_ic, get_motifs, calculate_motif_distance

test_meme_output = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/meme_bin/CLSTR1039_meme_out/'

records = get_motifs(test_meme_output, 1)

motif_a = records[2]
motif_b = records[3]

print(motif_a)

print(motif_b)

print(calculate_motif_distance(motif_a, motif_b))

'''
print(motif_a.instances[0])
print(motif_b.instances[0])


a_seqs = []
b_seqs = []


for offset in range(-len(motif_b) + 1, len(motif_a)):
    print('Working on offset ', offset, ':')
    get_ic(motif_a, motif_b, offset)'''
