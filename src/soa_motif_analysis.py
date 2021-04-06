'''
@author ichaudr

This module includes functions related to MEME and motif analysis, including executing MEME (Multiple Em for Motif Elicitation) calls through the command line and parses the output.
https://meme-suite.org/meme//doc/meme.html?man_type=web (command line version)

'''

from Bio.motifs import meme, Motif, Instances
from Bio import motifs
import os
import math
import itertools



def run_meme(input_file, output_dir, num_motifs=5, width_min=10, width_max=26, mode='anr', pal=False, meme_exec_path='/Users/ichaudr/Documents/UMBC/Lab-Erill/meme-5.3.0/src/meme'):
    '''
    Performs MEME analysis and saves output to a directory. The motifs can be parsed out by soa_meme_driver.get_motifs(). 

    Parameters
    ----------
    input_file: str
        Path to the FASTA file that has the sequences to be passed in for MEME analysis.
    output_dir: str
        Path to the output directory where all the MEME output will be saved. 
    pal: bool
        A MEME paramter that specifies if palindromes should be forced in the analysis. 
    meme_exec_path: str
        The path to the MEME executable

    num_motifs, width_min, width_max, mode -> see MEME command line documentation (https://meme-suite.org/meme//doc/meme.html?man_type=web)

    Returns
    ------
    None - output is stored in the output from the MEME command is stored in the specified output directory.

    '''

    #MEME command template
    meme_command = '{meme_exec_path} {input_file} -oc {output_dir} -nmotifs {num_motifs} -minw {width_min} -maxw {width_max} -mod {mode} -revcomp -dna'

    if pal:
        meme_command += " -pal"
    
    #Format and run MEME command
    os.system(meme_command.format(meme_exec_path=meme_exec_path, input_file=input_file, output_dir=output_dir, num_motifs=num_motifs, width_min=width_min, width_max=width_max, mode=mode))

def get_motifs(meme_data_dir, e_val_threshold):
    '''
    Parses the MEME output for the motifs.

    Parameters
    ----------
    meme_data_dir: str
        The path to the directory where the MEME results are stored (namely, the meme.xml file)
    e_val_threshold: float
        The maximum e-value for any motif that is returned. 
    
    Returns
    -------
    motifs: [Motifs]
        An array of Motif objects corresponding the motifs that were parsed from the MEME analysis.
    '''

    #List of Motif objects that met the threshold
    motifs_in_record = []

    #Pull all of the records from the MEME output file
    with open(meme_data_dir + 'meme.xml') as f:
        try:
            records = meme.read(f)
        except:
            print('Error with parsing MEME output.')
    
    if len(records) == 0:
        return motifs_in_record
    
    #Pull out motifs that meet the e value threhold
    for motif in records:
        if motif.evalue <= e_val_threshold:
            motifs_in_record.append(motif)
    
    return motifs_in_record 


################## Functions below are all collectively used to calculate the distance between two motifs. ###################
# Approach:
# 1. Determine the optimal alignment of the two motifs by maximizing the information content of the alignments.
# 2. Calculate the distance between each of the columns in the alignment.
# 3. The average distance between the columns is assigned to be the distance between the two motifs.

def get_alignment_offset(motif, other):
    '''
    Determines the optimal alignment of two motifs by maximizing the information content (ic) in the aligned regions.

    Parameters
    ----------
    motif, other: Motif objects
        The two motifs of interest.
    
    Returns
    -------
    offsets: int
        The offset that results in the maxium ic in the alignment. 
    '''
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
    '''
    Caculates the information content, ic, for a specific alignment. The approach makes a temporary motif object containing the overlapping sequences in the alignemnt and taking the average of the pssm.

    Parameters
    ----------
    motif, other: Motif objects
        The motifs of interest
    offset: int
        The offset value that results in the alignment of interest. 
    '''

    #Pull the sequences containined in the aligned region of the motifs from each of the motif instances. 
    alignment_len = min(len(motif)-offset, len(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in motif.instances]
    other_seqs = [site[:alignment_len] for site in other.instances]

    # Create the motif and compute the IC
    amotif = Motif(instances=Instances(motif_seqs+other_seqs))
    amotif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)

    #print('Motif Seqs: ' , motif_seqs)
    #print('Other Seqs: ' , other_seqs)
    #print('Offset ', offset)
    #print('IC: ' , amotif.pssm.mean(), '\n\n')

    return amotif.pssm.mean()

def calc_kld_distance(cola, colb):
    '''
    Calculates the KL distance between two columns of a pwm.

    Parameters
    ----------
    cola, colb: [float]
        Two columns A and B, repsectively, of a PWM. 

    Returns
    -------
    distance: float
        The calculated distance between the two columns. 
    '''

    safe_log2 = lambda x: math.log(x, 2) if x != 0 else 0.0
    distance = (sum(cola[l] * safe_log2(cola[l] / colb[l]) for l in "ACTG" if colb[l] != 0) + 
            sum(colb[l] * safe_log2(colb[l] / cola[l]) for l in "ACTG" if cola[l] != 0))
    return distance

def calc_euclidean(cola, colb):
    '''
    Calculates the euclidean distance between two columns of a pwm.

    Parameters
    ----------
    cola, colb: [float]
        Two columns A and B, repsectively, of a PWM. 

    Returns
    -------
    distance: float
        The calculated distance between the two columns. 
    '''
    distance = math.sqrt(sum((cola[let]-colb[let])**2 for let in "ACTG"))
    return distance

def pwm_col(motif_pwm, col):
    '''
    Extracts a specified column from a pwm. 

    Parameters
    ----------
    motif_pwm:
        The pwm from a motif object.
    
    col: int
        The column of interest.
    
    Returns
    ------
    col: dict {letter:frequency}
        A dictionary where the keys are A, T, C, and G and the values are their respective frequencies in the specified column. 
    '''
    col = dict((let, motif_pwm[let][col]) for let in "ACTG")
    return col

def calculate_motif_distance(motif, other, offset=None, padded=True, distance_function=calc_euclidean):
    '''
    Calculates the distance between two motifs by: (1) finding the maximum information content alignment and (2) determining the euclidian distance of that alignment.

    Parameters
    ----------
    motif, other: Motif objects
        The two motifs of interest.
    offset: int
        The alignment offset between the two motifs that results in the alignment with maximum information content.
    padded: bool
        If true, the distance calculation spans the entire length of both motifs and columns that are outside the alignment are compared to a column of equiprobable frequencies: ({A: 0.25, T: 0.25, C: 0.25, G: 0.25}). 
        If false, only the columns in the alignment are considered for the distance calculation.
    
    '''
    if offset is None:
        offset = get_alignment_offset(motif, other)
    if offset < 0:
        return calculate_motif_distance(other, motif, -1*offset, padded=padded, distance_function=distance_function)

    dists = []
    alignment_length = min(len(motif), len(other)-offset)

    for pos in range(alignment_length):
        cola = pwm_col(motif.pwm, pos)
        colb = pwm_col(other.pwm, pos+offset)

        #print(pos, '\tcolA: ', [seq[pos] for seq in motif.instances])
        #print(pos, '\tcolB', [seq[pos+offset] for seq in other.instances])

        dists.append(distance_function(cola, colb))
    
    if padded:
        #Add padded values for motif
        for pos in range(len(motif)):
            if pos > alignment_length - 1:
                cola = pwm_col(motif.pwm, pos)
                colb = {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}

                #print(pos, '\tpadded_1colA: ', [seq[pos] for seq in motif.instances])
                #print(pos, '\tpadded_1colB', ' ATCG')

                dists.append(distance_function(cola, colb))
        
        #Add padded values for other
        for pos in range(len(other) - 1):
            if pos < offset or pos > alignment_length + offset:
                cola = {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}
                colb = pwm_col(other.pwm, pos)

                #print(pos, '\tpadded_2colA: ', 'ATCG')
                #print(pos, '\tpadded_2colB', [seq[pos] for seq in other.instances])

                dists.append(distance_function(cola, colb))

    
    
    return sum(dists) / len(dists)


    ################## Functions below are all collectively used to determine whether a motif contains a direct or inverted repeat. ###################

def build_motif(sites):
    """Builds a Biopython motifs.Motif object out of given sites."""
    motif = motifs.create(sites)
    motif.pseudocounts = 0.8   # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
    return motif


def slice_sites(sites, start, end):
    """Slices each site."""
    return [site[start:end] for site in sites]


def score_site(pssm, site):
    """Scores the given site with the given PSSM."""
    return sum(pssm[site[i]][i] for i in range(len(site)))


def score_sites(pssm, sites):
    """Computes the average PSSM score of a list of sites."""
    return sum(score_site(pssm, site) for site in sites) / len(sites)


def direct_repeat(seq):
    """Match for a sequence in a direct-repeat pattern."""
    return seq


def inverted_repeat(seq):
    """Match for a sequence in an inverted-repeat pattern."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[b] for b in seq[::-1])


def find_pattern(motif, sites, self_score_ratio_threshold=0.6,
                 kmer_pair_score_ratio_threshold=0.3):
    """Finds pattern in a motif."""
    k = 4
    sites_len = len(sites[0])
    # Splices all sites into different sub-sequences (or k-mers) of size k
    all_kmers = [{'start': i, 'end': i+k, 'seqs': slice_sites(sites, i, i+k)}
                 for i in range(sites_len-k+1)]

    # Compute self-PSSM scores of the k-mers
    for kmer in all_kmers:
        motif = build_motif(kmer['seqs'])
        kmer['pssm'] = motif.pssm
        kmer['self_score'] = score_sites(kmer['pssm'], kmer['seqs'])

    # Find the highest self-PSSM score of all the k-mers
    max_self_score = max(kmer['self_score'] for kmer in all_kmers)
    # Keep only the k-mers with self_score values > self_score_ratio_threshold*max_self_score
    all_kmers = [kmer for kmer in all_kmers
                 if kmer['self_score'] > self_score_ratio_threshold*max_self_score]

    # Pattern is a tuple containing ('Type of repeat', score, start index, end index)
    # Set pattern to SB (single box) by default
    pattern = ('SB',0 , 0, 0)
    for kmer_a, kmer_b in itertools.combinations(all_kmers, 2):
        if not (kmer_a['start'] >= kmer_b['end'] or kmer_b['start'] >= kmer_a['end']):
            continue

        if kmer_a['self_score'] < kmer_b['self_score']:
            kmer_a, kmer_b = kmer_b, kmer_a

        # Look for direct-repeat
        score = score_sites(kmer_a['pssm'], kmer_b['seqs'])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and score > pattern[1]):
            pattern = ('DR', score, kmer_a['start'], kmer_b['start'])

        # Look for inverted-repeat
        score = score_sites(
            kmer_a['pssm'], [inverted_repeat(site) for site in kmer_b['seqs']])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and score > pattern[1]):
            pattern = ('IR', score, kmer_a['start'], kmer_b['start'])

    return pattern