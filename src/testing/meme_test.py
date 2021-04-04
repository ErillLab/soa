'''
Attempting to figure out the cause of this error:

EME analysis:  61%|█████████████████████████████████████████████████▉                                | 1285/2112 [2:37:23<1:41:17,  7.35s/it]
Traceback (most recent call last):
  File "sys_analysis.py", line 700, in <module>
    soa()
  File "sys_analysis.py", line 656, in soa
    cluster.motifs = get_motifs(meme_data_dir=out_dir, e_val_threshold=motif_e_val_threshold)
  File "/home/issac/Desktop/soa/src/soa_motif_analysis.py", line 71, in get_motifs
    records = meme.read(f)
  File "/home/issac/Desktop/soa/soa-venv/lib/python3.6/site-packages/Bio/motifs/meme.py", line 53, in read
    __read_motifs(record, xml_tree, sequence_id_name_map)
  File "/home/issac/Desktop/soa/soa-venv/lib/python3.6/site-packages/Bio/motifs/meme.py", line 176, in __read_motifs
    motif = Motif(record.alphabet, instances)
  File "/home/issac/Desktop/soa/soa-venv/lib/python3.6/site-packages/Bio/motifs/meme.py", line 67, in __init__
    motifs.Motif.__init__(self, alphabet, instances)
  File "/home/issac/Desktop/soa/soa-venv/lib/python3.6/site-packages/Bio/motifs/__init__.py", line 263, in __init__
    counts = self.instances.count()
  File "/home/issac/Desktop/soa/soa-venv/lib/python3.6/site-packages/Bio/motifs/__init__.py", line 219, in count
    counts[letter][position] += 1
KeyError: 'K'


This error was thrown when attempting to parse the results for CLSTR4. 
'''

from Bio.motifs import meme, Motif, Instances
from Bio import motifs

#Path to CLSTR4 meme output
meme_results = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/meme_bin/CLSTR4_meme_out/meme.xml'

with open(meme_results, 'r') as f:
    try:
      records = meme.read(f)
    except:
      print('error')
  
print('made it')
    