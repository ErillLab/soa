'''
@author ichaudr
A diagnostic script that was used to isolate BLAST searches that resulted incomplete BLAST hits. 
'''

import os
import pickle
from soa_features import AnnotatedHit

#Directory to the BLAST searches
blast_bin = '/Volumes/Issac_Ext/ErillsLab/SOA/blast_bin/'

blast_records = [file for file in os.listdir(blast_bin) if file.split('.')[-1]=='p']

for saved_record in blast_records:
    if not 'WP_001146752.1' in saved_record:
        pass

    curr_hits = pickle.load(open(blast_bin + saved_record, 'rb'))
    
    for hit in curr_hits:
        print(hit.five_end)

