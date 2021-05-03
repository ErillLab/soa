'''
A copy from https://github.com/ErillLab/Transfer_method_analysis/blob/71a042cc8ca0ce03a1d44c84d87f9bb509b6ec4c/src/motif.py#L81
Being used for testing purposes.
'''

from Bio.motifs import Motif, Instances
import math


test_motif = Motif(instances=Instances(["ATCAGTCA", "ATCAGTAA", "ATCTGTCA"]))
print(test_motif.consensus)