'''
@author ichaudr
A diagnostic script that was used to isolate BLAST searches that resulted incomplete BLAST hits. 
'''

import os
import pickle

#Directory to the BLAST searches
blast_bin = '/Volumes/Issac_Ext/ErillsLab/SOA/blast_bin'

blast_records = [file for file in os.listdir(blast_bin) if file.split('.')[-1]=='p']

print(blast_records)

