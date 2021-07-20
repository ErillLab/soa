'''
@author ichaudr

Holds functions relevant for similarity filter in sys_analysis.py.

'''

from Bio import pairwise2
from Bio.motifs import Motif, Instances
import random
import csv

def get_percent_matches(seq1, seq2, scoring={'match':2, 'mismatch':-1, 'gap_open':-1, 'gap_extd': -0.5}, cache_file=None):
    '''
    Returns the percent matches between two sequences. 

    Parameters
    ----------
    seq1, seq2: str
        The two sequences to calculate the percent matches for.
    scoring: dictionary
        The scoring to use for the alignment in the following format: scoring={'match':__, 'mismatch':__, 'gap_open':__, 'gap_extd': __}
        Default: scoring={'match':2, 'mismatch':-1, 'gap_open':-1, 'gap_extd': -0.5}
    cache_file:str
        The path to a temp file to store/pull alignment calculation results from. 
    
    Returns
    -------
    average_percent_matches: float
        The average percent match between seq1 and seq2. Calculated as: #matches / (#matches + #mismataches); gaps are not included
    '''

    if cache_file:
        with open(cache_file, 'r') as file:
            freader = csv.reader(file)
            for record in freader:
                if (record[0] == seq1 and record[1] == seq2) or (record[0] == seq2 and record[1] == seq1):
                    return record[2]
        

    #A list of percent matches of all the alignments - used to calculate the average 
    percent_matches = []
    for a in pairwise2.align.globalms(seq1, seq2, scoring['match'], scoring['mismatch'], scoring['gap_open'], scoring['gap_extd']):
        #The format of an alignment:
            # [0] : first sequence aligned
            # [1] : second sequenced aligned
            # [2] : alignment score
            # [3] : starting position of the alignment
            # [4] : ending position of the alignment 
    
        seq1_aligned = a[0]
        seq2_aligned = a[1]
        alignment_length = a[4]
        num_gaps = 0
        num_matches = 0

        for n in range(alignment_length):

            #Adjust the gap count
            if seq1_aligned[n] == '-' or seq2_aligned[n] == '-':
                num_gaps += 1
                continue

            #Adjust match count
            if seq1_aligned[n] == seq2_aligned[n]:
                num_matches += 1
        
        #Calculate the percent matches for the current alignment and append it to the list of percent matches
        percent_match = num_matches / (alignment_length - num_gaps)
        percent_matches.append(percent_match)
        
        #Return the average percent matches

    if len(percent_matches) <= 0:
        return 0

    average_percent_matches = sum(percent_matches)/len(percent_matches)

    if cache_file:
        with open(cache_file, 'a') as file:
            fwriter = csv.writer(file)
            fwriter.writerow([seq1, seq2,average_percent_matches])

    return average_percent_matches


def get_centroid(sequences):
    '''
    Determines, from a set of sequences, which sequence has the highest average percent identity to all of the other sequences in the set. (i.e. the centroid)

    Paramters
    ---------
    sequences: [str]
        A list of sequences
    
    Returns
    -------
    centroid: str
        The centroid sequences
    '''

    if len(sequences) == 0:
        return 0
    
    if len(sequences) <= 2:
        return sequences[0]

    #A dictionary holding the average percent id between a sequence and the rest of the sequences in the set {seq:average_percent_id}
    sim_list = {}

    for i in range(len(sequences)):
        
        #The percent similarities for this sequence
        percent_sims = []

        for j in range(len(sequences)):

            if i == j:
                continue

            percent_sims.append(get_percent_matches(sequences[i], sequences[j]))
        
        #Append the average percent similarity to the sim_list
        sim_list.update({sequences[i]:(sum(percent_sims)/len(percent_sims))})
    
    #Determine which sequence has the max average percent similarity and return it
    centroid = max(sim_list, key=sim_list.get)
    return centroid



def sim_filter(sequences, threshold_percent_id, cache_file=None):
    '''
    Filters a set of sequences so no two sequences have a percent ID higher than the threshold. 

    Parameters
    ----------
    sequences: [str]
        A list of sequences
    threshold_percent_id: float
        The threshold identity for the returned set of sequences (i.e. 0.75)
    cache_file:str
        The path to a temp file to store/pull alignment calculation results from. 

    Returns
    -------
    filtered_seqs: [str]
        A list of filtered sequences. 
    '''

    if threshold_percent_id == -1:
        return sequences

    if len(sequences) < 2:
        print('Not enough sequences passed in.')
        return sequences

    #Sort all of the sequences into "bins" based on whether or not their percent identity is above the threshold
    bins = {}

    i = 0
    while len(sequences) > 0:

        #Make a bin for the first sequence
        bin_key = "bin" + str(i)
        bins.update({bin_key:[sequences.pop(0)]})

        if len(sequences) == 0:
            continue

        #The list of sequnces that are within the threshold for the current bin
        matched_seqs =[]

        for seq in sequences:
            if len(seq) < 1:
                continue
            
            mp = float(get_percent_matches(seq1=seq, seq2=bins[bin_key][0], cache_file=cache_file))

            if mp == None:
                print(seq)
                continue

            if mp >= threshold_percent_id:
                matched_seqs.append(seq)
        
        #Append the matched seqs to the list for this bin and remove the seqs from the list of sequences
        for seq in matched_seqs:
            bins[bin_key].append(seq)
            sequences.remove(seq)
        
        i += 1
    
    #Get the centroid for each bin and return the final filtered list
    filtered_seqs = []
    for key in bins.keys():
        filtered_seqs.append(bins[key][0])
        #filtered_seqs.append(bins[key][int(random.random() * (len(bins[key]) - 1))])
        #filtered_seqs.append(get_centroid(bins[key]))
    
    return filtered_seqs