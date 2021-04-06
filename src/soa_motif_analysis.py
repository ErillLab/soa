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
import statistics
from soa_operon_cluster import OperonCluster


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

def calculate_motif_distance(motif, other, distance_function, offset=None, padded=True):
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
        return calculate_motif_distance(other, motif, distance_function, -1*offset, padded=padded)

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
# Approach:
# 1. Splice a given motif into kmers of a given length (i.e. 4).
# 2. Keep only the kmers with a high information content (filter out kmers that are noise).
# 3. Align kmers in all possible combinations and score each alignment.

def build_motif(sites):
    '''
    Builds a Biopython motifs.Motif object out of given sites.

    Parameters
    ----------
    sites
        A list of sequences that make up the motif
    
    Returns
    -------
    motif: motifs object
        The motif object created from sites 
    '''
    motif = motifs.create(sites)
    motif.pseudocounts = 0.8   # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
    return motif


def slice_sites(sites, start, end):
    '''
    Slices each site in sites[]

    Parameters
    ----------
    sites: [string]
        The list of sites (sequences) that make up a motif 
    start: int
        The index to begin splicing
    end: int
        The index to end splicing

    Returns
    -------
    spliced_sites: [string]
        The list of sites (sequences) spliced from start"end
    '''
    spliced_sites = [site[start:end] for site in sites]
    return spliced_sites


def score_site(pssm, site):
    '''
    Scores the given site with the given PSSM.

    Parameters
    ----------
    site: string
        sequence that is being scored with the given pssm (the pssm is not necessarily the site's pssm)
    pssm: [[float]]
        pssm that will be used to score site 

    Returns
    -------
    score: float
        The score of the site using the given pssm 
    '''
    score = sum(pssm[site[i]][i] for i in range(len(site)))
    return score


def score_sites(pssm, sites):
    '''
    Computes the average score of the sites[] passed in by score each site with the given pssm
    and then taking the average of all of the scores

    Parameters
    ----------
    sites: [string]
        sequences that are being scored with the given pssm (the pssm is not necessarily the sites' pssm)
    pssm: [[float]]
        pssm that will be used to score sites 

    Returns
    -------
    score: float
        The average score of each site in sites using the given pssm 
    '''
    average_score = sum(score_site(pssm, site) for site in sites) / len(sites)
    return average_score


def direct_repeat(seq):
    '''
    Returns the direct-repeat of the sequence passed in.

    Parameters
    ----------
    seq: string
        The seqence to be matched to a direct repeat pattern

    Returns
    -------
    seq: string
        The seqence matched to a direct repeat pattern
    '''
    return seq


def inverted_repeat(seq):
    '''
    Returns the inverted-repeat of the sequence passed in.

    Parameters
    ----------
    seq: string
        The seqence to be matched to a inverted repeat pattern

    Returns
    -------
    seq: string
        The seqence matched to an inverted repeat pattern
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[b] for b in seq[::-1])


def find_pattern(sites, self_score_ratio_threshold=0.6,
                 kmer_pair_score_ratio_threshold=0.3):
    '''
    Returns a boolean indicating whether a direct or inverted repeat pattern was found in a given motif

    Parameters
    ----------
    sites: [string]
        The sequences making up the motif
    self_score_ratio_threshold: float
        A threshold value to determine whether a score indicates a repeat
    kmer_pair_score_ratio_threshold: float
        A threshold value to keep only the kmers with a high information content

    Returns
    -------
    found_repeat: boolean
        Indicates whether a direct or inverted repeat was found
    '''
    # Splices all sites of a given motif into different sub-sequences (or k-mers) of size k
    k = 4
    sites_len = len(sites[0])
    all_kmers = [{'start': i, 'end': i+k, 'seqs': slice_sites(sites, i, i+k)}
                 for i in range(sites_len-k+1)]

    # Compute self-PSSM score of each k-mer
    for kmer in all_kmers:
        motif = build_motif(kmer['seqs'])
        kmer['pssm'] = motif.pssm
        kmer['self_score'] = score_sites(kmer['pssm'], kmer['seqs'])

    # Find the highest self-PSSM score of all the k-mers
    max_self_score = max(kmer['self_score'] for kmer in all_kmers)
    # Keep only the k-mers with self_score values > self_score_ratio_threshold * max_self_score
    all_kmers = [kmer for kmer in all_kmers
                 if kmer['self_score'] > self_score_ratio_threshold*max_self_score]

    # Pattern is a tuple containing ('Type of repeat', score, start index, end index)
    # Set pattern to SB (single box) and found_repeat to False by default
    pattern = ('SB',0 , 0, 0)
    found_repeat = False

    # Algin kmers in all possible combinations
    for kmer_a, kmer_b in itertools.combinations(all_kmers, 2):
        if not (kmer_a['start'] >= kmer_b['end'] or kmer_b['start'] >= kmer_a['end']):
            continue

        if kmer_a['self_score'] < kmer_b['self_score']:
            kmer_a, kmer_b = kmer_b, kmer_a

        # Look for a direct 
        # Score kmerB using the pssm of kmerA, if the score > threshold, a repeat has been found
        score = score_sites(kmer_a['pssm'], kmer_b['seqs'])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and score > pattern[1]):
            pattern = ('DR', score, kmer_a['start'], kmer_b['start'])
            found_repeat = True

        # Look for inverted-repeat
        score = score_sites(
            kmer_a['pssm'], [inverted_repeat(site) for site in kmer_b['seqs']])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and score > pattern[1]):
            pattern = ('IR', score, kmer_a['start'], kmer_b['start'])
            found_repeat = True

    return found_repeat

################## Functions below are used to create csv files containing edges and nodes for a Gephi network. ###################
# Approach:
# 1. Filter motifs returned by MEME, keeping only those with direct and inverted repeats.
# 2. Calculate the KLD distances between each motif in two given operons, and report the distance between two operons as the:
#       a. average of all the distances between all motifs in both operons
#       b. minimum of all the distances between all motifs in both operons
#       c. maximum of all the distances between all motifs in both operons
#       d. median of all the distances between all motifs in both operons
#       e. noise reduction distance calculation
#           The operon with fewer motifs is Operon A, and the other is Operon B. For each motif in Operon A, 
#           find the motif in Operon B for which the motif in Operon A has the smallest distance. Then, take
#           the average of these minimum distances to calculate the distance between the operons
# 3. Write to nodes and edges csv files for each different distance calculation

def get_operon_clusters():
    '''
    Parses the meme_bin folder and creates operon_cluster objects for each cluster directory

    Parameters
    ----------
    None
    
    Returns
    -------
    operon_clusters:[operon_cluster] 
        An array holding operon_cluster objects for each cluster in meme_bin
    '''

    # An array of operon cluster objects
    operon_clusters = []

    # Get the file path for meme_bin
    local_path = os.path.dirname(os.getcwd())
    clusters_path = local_path + '/meme_bin/'

    # The maximum e-value for any motif that is returned by get_motifs()
    motif_e_val_threshold = 10E-3

    # Traverse the meme_bin directory, create operon_cluster objects and set their motifs after filtering for repeats
    for cluster_dir in os.listdir(clusters_path):
        if(cluster_dir[0:5] == "CLSTR"):
            print(cluster_dir)
            # Create operon cluster object, get operon ID from file name
            new_cluster = OperonCluster(cluster_id = cluster_dir.split("_")[0][5:])
            # Get all motifs associated with the cluster
            cluster_motifs =  get_motifs(meme_data_dir = clusters_path + cluster_dir + '/', e_val_threshold=motif_e_val_threshold)

            # Keep only the motifs with direct or inverted repeats
            motifs_with_repeats = []
            for motif in cluster_motifs:
                if (find_pattern(motif.instances)):
                    motifs_with_repeats.append(motif)
            new_cluster.motifs = motifs_with_repeats
            
            # Add new_cluster object to the operon_clusters[] array
            operon_clusters.append(new_cluster)

    return operon_clusters

def write_to_edge_file(file_handle, cluster_one, cluster_two, value):
    '''
    Writes cluster_one, cluster_two, edge_value to edge csv file

    Parameters
    ----------
    file_handle: fileObj
        The file that is being written to

    cluster_one, cluster_two: string
        The ID's of the two clusters who's distances are being calculated
    
    value: float
        The value being written to the csv file for the edge (average, median, min, or max)
    
    Returns
    -------
    None - the file is written to by the file handles
    '''
    file_handle.write(cluster_one)
    file_handle.write(',')
    file_handle.write(cluster_two)
    file_handle.write(',')
    file_handle.write(value)
    file_handle.write('\n')

def write_to_csv(filename, operon_clusters, distance_calc):
    '''
    Creates and writes edge and node csv files

    Parameters
    ----------
    operon_clusters:[operon_cluster] 
        An array holding operon_cluster objects for each cluster in meme_bin

    Returns
    -------
    None - writes distance and node information to csv files
    '''

    # Create file objects for various distance calculation methods (for edge csv files)
    edges_output_file_average = open(filename + "_edges_average.csv", "w")
    edges_output_file_median = open(filename + "_edges_median.csv", "w")
    edges_output_file_minimum = open(filename + "_edges_minimum.csv", "w")
    edges_output_file_maximum = open(filename + "_edges_maximum.csv", "w")
    edges_output_file_noise_reduction = open(filename + "_edges_noise_reduction.csv", "w")
    
    #  Write headers to be recognized in Gephi
    edges_output_file_average.write("Source,Target,Weight\n")
    edges_output_file_median.write("Source,Target,Weight\n")
    edges_output_file_minimum.write("Source,Target,Weight\n")
    edges_output_file_maximum.write("Source,Target,Weight\n")
    edges_output_file_noise_reduction.write("Source,Target,Weight\n")

    # Create file oject for node csv file
    nodes_output_file = open(filename + "_nodes.csv", "w")
    nodes_output_file.write("Id,Label\n")

    for i in range(0,len(operon_clusters)):
        # Create a node for every cluster in nodes.csv
        nodes_output_file.write(operon_clusters[i].cluster_id)
        nodes_output_file.write(',')
        nodes_output_file.write(operon_clusters[i].cluster_id)
        nodes_output_file.write("\n")

        # Compare the motifs of each cluster with the motifs of the other cluster
        for j in range(i,len(operon_clusters)):
            # Don't compare a cluster with itself
            if (i != j):
                # Hold all pairwise distances and minimum distances
                pairwise_distances = []
                minimum_distances = []
                # Identify the operon with fewer motifs
                # Always search for the minimum distance between a motif from the operon with the 
                # fewer amount of motifs and all of the motifs of the operon with the larger amount
                # of motifs (noise reduction)
                if(len(operon_clusters[i].motifs) <= len(operon_clusters[j].motifs)):
                    operon_fewer_motifs = operon_clusters[i].motifs
                    operon_more_motifs = operon_clusters[j].motifs
                else:
                    operon_fewer_motifs = operon_clusters[j].motifs
                    operon_more_motifs = operon_clusters[i].motifs

                # Compare all of the motifs in cluster one with the motifs of cluster two
                for k in range(0, len(operon_fewer_motifs)):
                    # Set initial minimum value 
                    if (len(operon_more_motifs) > 0):
                        min_val = calculate_motif_distance(operon_fewer_motifs[k],operon_more_motifs[0], distance_calc)
                    
                    for l in range(0, len(operon_more_motifs)):
                        distance_between_motifs = calculate_motif_distance(operon_fewer_motifs[k],operon_more_motifs[l], distance_calc)
                        #print("fewer: (",len(operon_fewer_motifs),") ",k, " more: (",len(operon_more_motifs),") ", l, " distance: ", distance_between_motifs)
                        pairwise_distances.append(distance_between_motifs)

                        # If found new minimum, set min_val as new minimum value
                        if (distance_between_motifs < min_val):
                            min_val = distance_between_motifs

                        # If last comparison, add min_val to minimum_distances array
                        if (l == len(operon_more_motifs)-1):
                            minimum_distances.append(min_val)
                            #print("min_val: ", min_val)
                        
                # If minimum distances for noise reduction were calculated
                if (len(minimum_distances) > 0):
                    write_to_edge_file(edges_output_file_noise_reduction,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(sum(minimum_distances)/len(minimum_distances)))

                # if pairwise_distances were calculated, write average, median, min, and max values to csv filie
                if len(pairwise_distances) > 0:
                    write_to_edge_file(edges_output_file_average,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(sum(pairwise_distances)/len(pairwise_distances)))
                    write_to_edge_file(edges_output_file_median,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(statistics.median(pairwise_distances)))
                    write_to_edge_file(edges_output_file_minimum,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(min(pairwise_distances)))
                    write_to_edge_file(edges_output_file_maximum,operon_clusters[i].cluster_id, operon_clusters[j].cluster_id, str(max(pairwise_distances)))
                    print(i+1, "/" , len(operon_clusters), "--", operon_clusters[i].cluster_id, "|",operon_clusters[j].cluster_id)

    # close all output files
    edges_output_file_average.close()
    edges_output_file_median.close()
    edges_output_file_maximum.close()
    edges_output_file_minimum.close()
    edges_output_file_noise_reduction.close()
    nodes_output_file.close()

if __name__ == "__main__":
    distance_calc = calc_euclidean
    all_clusters = get_operon_clusters()
    write_to_csv(all_clusters, distance_function=distance_calc)