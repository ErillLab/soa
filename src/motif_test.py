'''
A copy from https://github.com/ErillLab/Transfer_method_analysis/blob/71a042cc8ca0ce03a1d44c84d87f9bb509b6ec4c/src/motif.py#L81
Being used for testing purposes.
'''

from Bio import motifs
import math

def alignment(motif, other):
    max_ic = float('-inf')
    for offset in range(-len(motif) + 1, len(other)):
        if offset < 0:
            ic = ic_at(motif, other, -offset)
        else:
            ic = ic_at(other, motif, offset)

        
        if ic > max_ic:
            max_ic = ic
            max_offset = offset
    return max_offset

def ic_at(motif, other, offset):
    alignment_len = min(len(motif)-offset, len(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in motif.instances]
    other_seqs = [site[:alignment_len] for site in other.instances]

    # Create the motif and compute the IC
    amotif = motifs.Motif(instances=motifs.Instances(motif_seqs+other_seqs))
    amotif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)


    print('Seqs from motif for offset=', offset, ':  ', motif_seqs)
    print('Seqs from other for offset=', offset, ':  ', other_seqs)
    print(amotif.pssm.mean())
    print('---')

    return amotif.pssm.mean()

def euclidean(cola, colb):
    """Euclidean distance between two frequency vector"""
    return math.sqrt(sum((cola[let]-colb[let])**2 for let in "ACTG"))

def pwm_col(motif, col):
    """A column of the PWM"""
    return dict((let, motif.pwm[let][col]) for let in "ACTG")

def distance(motif, other, col_dist_func=euclidean, offset=None):
    """Compute the distance between two motifs, using the given column distance
    function."""
    if offset is None:
        offset = alignment(motif, other)
    if offset < 0:
        return distance(other, motif, col_dist_func, -offset)
    dists = []
    print('offset: ', offset)
    for pos in range(min(len(motif), len(other)-offset)):
        cola = pwm_col(motif, pos)
        colb = pwm_col(other, pos+offset)
        dists.append(col_dist_func(cola, colb))
    return sum(dists) / len(dists)