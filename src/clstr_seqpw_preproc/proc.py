'''
@author ichaudr

The goal of the script is to compute percent matches for all sequence pairs for each cluster outputted from the SOA pipeline. This will speed up the process for filtering the promoter sequences. 

'''

import os
import sys
sys.path.append(os.path.abspath('../'))
from soa_sim_filter import get_percent_matches
from Bio import SeqIO
from tqdm import tqdm
import json

final_data = {}

#Path to all promoters
prom_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/CLSTR_Promoters_No_Filter_05262021/'
promoter_files = [f for f in list(os.listdir(prom_path)) if ".fasta" in f]

for f in promoter_files[0:2]:
    full_file = prom_path + f
    print(full_file)

    #Read in records
    records = list(SeqIO.parse(full_file, 'fasta'))[0:5]

    #For each record, determine all sequence-pair percent matches and append to final_data
    for r1 in records:
        curr_matches = {}
        for r2 in tqdm(records, desc=r1.id):
            curr_match = get_percent_matches(seq1=str(r1.seq), seq2=str(r2.seq))
            curr_matches[str(r2.id)] = curr_match
        final_data[str(r1.id)] = curr_matches
    print('-'*50)

json.dump(final_data, fp=open('./processed_seq_matches.json', 'w'))