'''
@author ichaudr

'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from soa_operon import Operon
from soa_sim_filter import sim_filter
from tqdm import tqdm


class OperonCluster:

    def __init__(self, cluster_id):
        self.cluster_id = cluster_id
        self.operons = []
        self.filtered_promoters = []
        self.motifs = []
    
    def filter_promoters(self, threhold_percent_id=0.85):
        '''
        Pulls the promoter from all of the operons associated with this cluster and filters them based off the similarity threshold. 

        Parameters
        ----------
        threshold_similarity: float
            The maximum similarity allowed between any two promoters.

        Returns
        -------
        None
        '''

        #Pool all of the promoters from the operons in this cluster
        all_promoters = []

        for op in self.operons:
            all_promoters.append(op.promoter)

        tqdm.write('----- Pre-filtering: ' + str(len(all_promoters)))
        
        #Filter the promoters
        self.filtered_promoters = sim_filter(all_promoters, threhold_percent_id=threhold_percent_id)

        tqdm.write('----- Post-filtering: ' + str(len(self.filtered_promoters)))


    def write_promoters(self, output_file):
        '''
        Writes the filtered promoters to an output FASTA file to be used in the MEME motif discovery.

        Parameters
        ---------
        output_file: str
            File path and name to write to
        
        Returns
        ------
        None 
        '''

        to_write = []

        #Create SeqRecord objects for all of the filtered promoters
        i = 0
        for prom in self.filtered_promoters:
            seq_id = self.cluster_id + '_p' + str(i)
            to_write.append(SeqRecord(Seq(prom),id=seq_id, description='|'))
            i += 1
        
        #Write all sequences to file
        SeqIO.write(to_write, output_file, 'fasta')

    def __str__(self):
        '''
        Overriden str()
        '''

        to_return = 'Operon Cluster: ' + str(self.cluster_id) + '\nNumber of operons: ' + str(len(self.operons)) + '\nMotifs:\n'

        for m in self.motifs:
            to_return += str(m)
            to_return += '\n'
        
        return to_return