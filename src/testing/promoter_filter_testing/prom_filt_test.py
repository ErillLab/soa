'''
@author ichaudr
May 26, 2021

The purpose of this script is to analyze the output of the sequence filtering step at different similarity thresholds. 
'''

from Bio import SeqIO
import os
import sys

import numpy as np
import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from scipy.stats import norm
import statistics as st
from tqdm import tqdm

from Bio import pairwise2
import random
import csv



def plthisto_and_draw(data=[], file='', title='', resolution = 0.01):
    '''
    A function to draw a histogram distribution of the data and save the plot to file. 
    '''


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(data, bins=len(data))
    ax.set_ylim([0, 20])

    title += '\nMean: ' + str(st.mean(data)) + '\nStd Dev: ' + str(st.mean(data))+ '\nMax: ' + str(max(data)) + '\nMin: ' + str(min(data))

    ax.set_title(title, fontsize=8)
    ax.set_xlabel('Number of Promoter Sequences', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)

    plt.savefig('histo_'+ file)


def pltnorm_and_draw(data=[], file='', title='', resolution = 0.01):
    '''
    A function to draw a normal distribution of the data and save the plot to file. 
    '''

    x_axis = np.arange(0, max(data), resolution)
    mean = st.mean(data)
    sd = st.stdev(data)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_axis, norm.pdf(x_axis, mean, sd))

    title += '\nMean: ' + str(mean) + '\nStd Dev: ' + str(sd) + '\nMax: ' + str(max(data)) + '\nMin: ' + str(min(data))

    ax.set_title(title, fontsize=8)
    ax.set_xlabel('Number of Promoter Sequences', fontsize=12)

    plt.savefig('norm_' + file)


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

def sim_filter(sequences, threshold_percent_id, cache_file='default_cache.csv'):
    '''
    Filters a set of sequences so no two sequences have a percent ID higher than the threshold. 

    Parameters
    ----------
    sequences: [str]
        A list of sequences
    threshold_percent_id: float
        The threshold identity for the returned set of sequences (i.e. 0.75)

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

            mp = float(get_percent_matches(seq, bins[bin_key][0], cache_file=cache_file))

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
        filtered_seqs.append(bins[key][int(random.random() * (len(bins[key]) - 1))])
        #filtered_seqs.append(get_centroid(bins[key]))
    
    return filtered_seqs


#Path to all promoters
prom_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/CLSTR_Promoters_No_Filter_05262021/'
promoter_files = [prom_path+f for f in list(os.listdir(prom_path)) if ".fasta" in f]

filter_levels = np.arange(.50, 1.01, .10)

#Plot number of promoters without filtering
prom_counts = []
for f in tqdm(promoter_files, desc='Making Plot With No Filtering'):
    records = list(SeqIO.parse(f,'fasta'))
    seqs = [str(r.seq) for r in records]
    to_add = len(seqs)
    prom_counts.append(to_add)


pltnorm_and_draw(prom_counts, title='No Filter', file='no_filter.png')
plthisto_and_draw(prom_counts, title='No Filter', file='no_filter.png')


#Plot number of promoters with filtering
for t in filter_levels:
    prom_counts = []
    for f in tqdm(promoter_files, desc='Filtering at '+ str(t)):
        cache_file = str(f.split('/')[-1]) + '_prom_sim_cache.csv'
        open(cache_file, 'a')
        records = list(SeqIO.parse(f,'fasta'))
        seqs = [str(r.seq) for r in records]
        tqdm.write('\t' + f + '\t' + str(len(seqs)))
        to_add = len(sim_filter(sequences=seqs, threshold_percent_id=t, cache_file=cache_file))
        prom_counts.append(to_add)

    pltnorm_and_draw(prom_counts, title=('Threshold Similarity Level: ' + str(t)), file=('filter_' + str(t) + '.png'))
    plthisto_and_draw(prom_counts, title=('Threshold Similarity Level: ' + str(t)), file=('filter_' + str(t) + '.png'))  





