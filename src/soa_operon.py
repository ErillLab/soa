'''
@author: ichaudr

'''

from soa_features import AnnotatedHit, GenomeFeature

class Operon:

    def __init__(self, operon_id, genome_fragment_name, genome_accession, genome_features, strand):
        self.features = []
        self.operon_id = operon_id
        self.cluster_id = operon_id.split('_')[3].replace(' ', '')
        self.genome_fragment_name = genome_fragment_name
        self.genome_accession = genome_accession
        self.genome_features = genome_features
        self.strand = strand
        self.promoter = None
    
    def add_feature(self, feature):
        '''
        Appends feature to the list of features associated with this operon and sorts it from 5' to 3'

        Parameters
        ----------
        feature: GenomeFeature object
            Feature to be added
        
        Returns
        -------
        None
        '''
        #Checking if the feature is present in the features associated with the genome for this operon.
        if feature in self.genome_features:
            self.features.append(feature)
            self.features = sorted(self.features, key=lambda feature: feature.five_end)
        else:
            raise Exception("The feature you are trying to add is not in the genome assigned for this operon.")
    
    
    
    def __str__(self):
        try:
            to_return = "OPERON:" + str(self.genome_fragment_name) + "(" + str(self.genome_accession) + ")\nStrand: " + str(self.strand) + "\nOperon_id: " + str(self.operon_id) + "\n"

            for f in self.features:
                if isinstance(f, AnnotatedHit):
                    to_return = to_return + "\tOriginal Query: " + f.query_accession + "\tHit Accession: " + f.protein_accession + "\tLocus Tag: " + f.locus_tag + "\n"
                elif isinstance(f, GenomeFeature):
                    to_return = to_return + "\tIntergenic Feature:  " + f.protein_accession + "\tLocus Tag: " + f.locus_tag + "\n"
        except:
            to_return = "An error occured printing an operon for " + str(self.genome_fragment_name) +  "(" + str(self.genome_accession) + ")\n"
        
        return to_return