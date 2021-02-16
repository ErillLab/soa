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

