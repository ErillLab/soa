'''
@author ichaudr

This module includes functions related to MEME and motif analysis, including executing MEME (Multiple Em for Motif Elicitation) calls through the command line and parses the output.
https://meme-suite.org/meme//doc/meme.html?man_type=web (command line version)

'''

from Bio.motifs import meme
import os



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
    motifs = []

    #Pull all of the records from the MEME output file
    with open(meme_data_dir + 'meme.xml') as f:
        records = meme.read(f)
    
    #Pull out motifs that meet the e value threhold
    for motif in records:
        if motif.evalue <= e_val_threshold:
            motifs.append(motif)
    
    return motifs 



################## Functions below are all collectively used to calculate the distance between two motifs. ###################
# Approach:
# 1. Determine the optimal alignment of the two motifs by maximizing the information content of the alignments.
# 2. Calculate the distance between each of the columns in the alignment.
# 3. The average distance between the columns is assigned to be the distance between the two motifs.

def get_ic(motif_a, motif_b, offset):
    '''
    Calculates the information content for the alignment of interest. 

    Parameters
    ---------
    motif_a, motif_b: Motif objects
        The two motifs of interest.
    offset: int
        The alignment is represented as an offset value. See soa_motif_analysis.get_alignment_offset(). 
    '''

    #Pull out the sequences from motif_a that are part of the alignemt
    a_seqs = []

    for sequence in motif_a.instances:

        start_pos = 0
        end_pos = len(sequence)

        if offset < 0:
            start_pos = 0
        
        if offset > 0:
            start_pos = offset
        
        end_pos = len(motif_b) + offset
        a_seqs.append(sequence[start_pos:end_pos])
    
    print(a_seqs[0])
    

    #Pull out the sequences from motif_b that are part of the alignment
    b_seqs = []

    for sequence in motif_b.instances:
        
        start_pos = 0
        end_pos = len(sequence)

        if offset < 0:
            
            if len(motif_a) == len(motif_b):
                start_pos = -(len(motif_b) + offset)
                end_pos = None
            else:
                start_pos = offset * -1
                end_pos = len(motif_a) + offset
        
        if offset > 0:
            start_pos = 0
            end_pos = len(motif_a) - offset
        
        
        b_seqs.append(sequence[start_pos:end_pos])

    print(b_seqs[0])



def get_alignment_offset(motif_a, motif_b):
    '''
    Determines the best alignment between the two motifs and returns the offset required to obtain the alignment. The optimal alignment is chosen as the gapless alignment that maximizes the information content shared between the two motifs. 

    Parameters
    ----------
    motif_a, motif_b: Motif objects
        The two motifs of interest.
    
    Returns
    -------
    alignment_offset: [int]
        The alignment is returned as a list equivalent offeset values. A offset value is the number of columns between the first position of motif_a and motif_b. The offset has range: -len(motif_a) + 1 < offset < len(motif_b).
        At the minimum value of -len(motif_a) + 1, the offset means that the first position of motif_b is aligned with the last position in motif_a. 
        At the maximum value of len(motif_b), the offset means that the first position of motif_a is aligned with the last position of motif_a.
        An offset of zero means that the first position of motif_a is aligned with the first position of motif_b.
    '''

    #Holds the maximum information content. Currently set to negative infinity. 
    max_ic = float('-inf')
 
    alignment_offsets = []

    for offset in range(-len(motif_b) + 1, len(motif_a)):
        
        #Calculate the information content for the alignment with the current offset. 
        curr_ic = 0

        if curr_ic >= max_ic:
            alignment_offsets.append(offset)
    
    return alignment_offsets



def calculate_motif_distance(motif_one, motif_two):
    '''
    Returns the distance between two motifs by:
        1. Finding the optimal alignment that maximizes the information content. 
        2. Computes the Euclidian distance or the KL divergence between the aligned columns. 
    
    Parameters
    ----------
    motif_one, motif_two: Motif objects

    Returns
    -------
    motif_distance: float
        The computed distance between the two motifs.
    '''

    pass

