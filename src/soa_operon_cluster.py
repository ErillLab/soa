'''
@author ichaudr

'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.motifs import Motif, Instances
from soa_operon import Operon
from soa_sim_filter import sim_filter
from tqdm import tqdm
import json


class OperonCluster:

    def __init__(self, cluster_id=None):
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

    def filter_motifs_by_instance_count(self, instances_min=5):
        '''
        Filters out any motifs that do not meet the threhold number of instances.

        Parameters
        ----------
        instances_min: int
            The minimum number of instances required for a motif to remain in the cluster. 
        '''
        self.motifs = [m for m in self.motifs if len(m.instances) >= instances_min]

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
    
    def export_to_json(self, output_file):
        '''
        Exports all the information for this cluster into a json file.

        Parameters
        ----------
        output_file: str
            The JSON file to write this cluster to.
        
        Returns
        -------
        None
        '''

        #A dictionary holding all of the data for this cluster
        data = {
            "cluster_id":self.cluster_id,

            "operons":[{

                "operon_id":op.operon_id,
                "cluster_id":op.cluster_id,
                "genome_accession":op.genome_accession,
                "genome_fragment_name":op.geonome_fragment_name,
                "features":[feat.protein_accession for feat in op.features],
                "promoter":op.promoter,
                "strand":op.strand

            } for op in self.operons],

            "filtered_promoters":self.filtered_promoters,
            "motifs":[str(m.instances).split('\n')[:-1] for m in self.motifs]
        }

        with open(output_file, 'w') as file:
            json.dump(data, file)


    def load_from_json(self, file_path):
        '''
        Sets up all memeber variables based on the stored data in the JSON file. 

        Parameters
        ----------
        file_path: str
            Path to the JSON file that will be loaded in.
        
        Returns
        -------
        None
        '''

        file_reader = json.load(open(file_path, 'r'))

        self.cluster_id = file_reader['cluster_id']

        for op in file_reader['operons']:
            temp_op = Operon(
                operon_id=op['operon_id'],
                genome_fragment_name='imported_cluster',
                genome_accession=op['genome_accession'],
                genome_features='imported',
                strand='/'
            )

            temp_op.features=op['features']

            temp_op.promoter = op['promoter']

            self.operons.append(temp_op)
        
        self.filtered_promoters = file_reader['filtered_promoters']
        
        for m in file_reader['motifs']:
            self.motifs.append(Motif(instances=Instances(m)))


    def __str__(self):
        '''
        Overriden str()
        '''

        to_return = 'Operon Cluster: ' + str(self.cluster_id) + '\nNumber of operons: ' + str(len(self.operons)) + '\nMotifs:\n'

        for m in self.motifs:
            to_return += str(m)
            to_return += '\n'
        
        return to_return