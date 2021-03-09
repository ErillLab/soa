'''
@author: ichaudr

'''

from Bio import Entrez
import time


class GenomeFeature:
    '''
    Holds information for a feature in a gene record. Specifically carries the following information:
        - Genome accession number
        - Coding start and coding stop positions for the feature
        - The strand the gene is located on, + or -
        - The protein ID for the protein product of that feature
    '''

    def __init__(self, genome_accession, genome_fragment_name, req_limit, sleep_time, strand, aa_sequence=None, coding_start=None, coding_end=None, five_end=-1, three_end=-1, protein_accession="None", locus_tag="None"):
        self.genome_accession = genome_accession
        self.protein_accession = protein_accession
        self.locus_tag = locus_tag
        self.genome_fragment_name = genome_fragment_name
        self.req_limit = req_limit
        self.sleep_time = sleep_time
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.strand = strand
        self.five_end = five_end
        self.three_end = three_end
        self.aa_sequence = aa_sequence
    
    def get_intergenic_distance(self, other):
        '''
        Determines the distance to another feature object. 

        Parameters
        ----------
        other: GenomeFeature object
            The feature to be compared to this feature.
        
        Return
        ------
        distance: int
            The absolute value of the distance between the two features.
        '''


        distance = int(min(abs((int(self.three_end) - int(other.five_end))), abs((int(self.five_end) - int(other.three_end)))))
        return distance


    def __str__(self):
        to_return = ""
        to_return = to_return + "Genome Feature\n"
        to_return = to_return + str(self.genome_fragment_name) + "\n"
        to_return = to_return + str(self.genome_accession) + ": Coding start position: " + str(self.coding_start) + " Coding end position: " + str(self.coding_end) + "\n"
        to_return = to_return + "Strand: " + str(self.strand) + "\n"
        to_return = to_return + "Protein ID: " + str(self.protein_accession) + "\n"
        to_return = to_return + "Locus Tag: " + str(self.locus_tag) + "\n"

        return to_return
    
    def __eq__(self, other):
        if other == None:
            return False

        if (self.coding_start == other.coding_start and self.coding_end == other.coding_end and self.strand == other.strand) or self.locus_tag == other.locus_tag:
            return True
        else:
            return False


class AnnotatedHit(GenomeFeature):
    '''
    Holds the annotated information for a BLAST hit. It an extended form of a feature that also holds infromation
    about the alignment start/stop positions and the query accesion. 
    
    '''
    def __init__(self,query_accession, operon_id, hit_accession, genome_fragment_name, align_start, align_end, strand, alignment_seq, percent_identity, req_limit, sleep_time):
        self.query_accession = query_accession
        self.operon_id = operon_id
        self.align_start = align_start
        self.align_end = align_end
        self.feature_found = False
        self.percent_identity = percent_identity
        self.alignment_seq = alignment_seq
        
        if strand > 0:
            self.strand = '+'
        elif strand < 0:
            self.strand = '-'

        super().__init__(genome_accession=hit_accession, genome_fragment_name=genome_fragment_name, req_limit=req_limit, sleep_time=sleep_time, strand=self.strand)
        
        
    def fetch_feature(self, record, margin_limit=20, max_attempts=5, mult_factor=3):
        '''
        Determines the exact feature that corresponds to the alignment range in the genome. It will pull the entire annotated genome of the genome fragmet and search through
        all the features to determine which one contains the alignmet. It will than recalculate the percent identitiy by aligning the entire feature with the reference sequence. 

        Parameters
        ----------
        record: XML parsed object
            The full genome record that this hit belongs to. 
        margin_limit: int
            The max "wiggle" room between the 5' and 3' bounds of a feature and the bounds of the alignment for that feature to be assigned to this hit.
        max_attemps: int
            The max number of times to attempt finding the feature
        mult_factor: int
            The margin_limit will be multiplied by this parameter every through every attempt.

        Returns
        -------
        None - the information is stored in the object

        '''
        '''
        #####TESTING#####

        try:
            handle = Entrez.efetch(db="nuccore", id= self.genome_accession, strand=self.strand, seq_start=self.align_start, seq_stop=self.align_end, rettype='gbwithparts', retmode='XML')
            record_1 = Entrez.read(handle, 'xml')
        except:
            print('error')

        locus_id = ''
        protein_id = ''

        for feature in record_1[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'CDS':
                for quality in feature['GBFeature_quals']:
                    #Check for protein id
                    if quality['GBQualifier_name'] == 'protein_id':
                        protein_id = quality['GBQualifier_value']
                    
                    #Check for locus tag
                    if quality['GBQualifier_name'] == 'locus_tag':
                        locus_id = quality['GBQualifier_value']
        
        for feature in record[0]['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'CDS':
                    t_pid = ''
                    t_lid = ''

                    for quality in feature['GBFeature_quals']:
                        if quality['GBQualifier_name'] == 'protein_id':
                            t_pid = quality['GBQualifier_value']
                        if quality['GBQualifier_name'] == 'locus_tag':
                            t_lid = quality['GBQualifier_value']
                            
                    
                    if t_lid == locus_id and t_pid == protein_id:
                        self.locus_tag = t_lid
                        self.protein_accession = t_pid
                        self.feature_found = True
                        self.coding_start = int(feature['GBFeature_intervals'][0]['GBInterval_from'])
                        self.coding_end = int(feature['GBFeature_intervals'][0]['GBInterval_to'])
                        self.five_end = min(self.coding_end, self.coding_start)
                        self.three_end = max(self.coding_start, self.coding_end)



        #################
        '''        
        #Set the 5' and 3' bounds of the alignment
        align_five_end = min(self.align_start, self.align_end)
        align_three_end = max(self.align_start, self.align_end)

        continue_searching = True
        current_num_attempts = 1

        
        while continue_searching:

            if current_num_attempts > max_attempts:
                continue_searching = False
                continue

            #Pulling all features from the Entrez results
            for feature in record[0]['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'CDS':
                    if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                        coding_start = min(int(feature['GBFeature_intervals'][0]['GBInterval_from']), int(feature['GBFeature_intervals'][0]['GBInterval_to']))
                        coding_end = max(int(feature['GBFeature_intervals'][0]['GBInterval_from']), int(feature['GBFeature_intervals'][0]['GBInterval_to']))

                        #Determining if feature is inclusive of the alignment. 
                        if (abs(min(coding_start, coding_end) - align_five_end) < margin_limit) and (abs(max(coding_start, coding_end) - align_three_end) < margin_limit):

                            self.feature_found = True

                            self.coding_start = coding_start
                            self.coding_end = coding_end
                            
                            #Parse the protein ID, locus tag, and protein sequence from the feature table
                            for quality in feature['GBFeature_quals']:
                                if quality['GBQualifier_name'] == 'protein_id':
                                    protein_accession = quality['GBQualifier_value']
                                    if protein_accession == None:
                                        self.protein_accession = "None"
                                    else:
                                        self.protein_accession = protein_accession

                                if quality['GBQualifier_name'] == 'locus_tag':
                                    locus_tag = quality['GBQualifier_value']
                                    if locus_tag == None:
                                        self.locus_tag = "None"
                                    else:
                                        self.locus_tag = locus_tag
                                
                                if quality['GBQualifier_name'] == 'translation':
                                    sequence = quality['GBQualifier_value']
                                    if sequence == None:
                                        self.aa_sequence = 'None'
                                    else:
                                        self.aa_sequence = sequence

                            self.five_end = min(coding_end, coding_start) 
                            self.three_end = max(coding_start, coding_end)

                            continue_searching = False
                            current_num_attempts = current_num_attempts + 1
            
            if not self.feature_found:
                current_num_attempts = current_num_attempts + 1
                margin_limit = margin_limit * mult_factor
                

        if not self.feature_found:
            print("Error: No feature found for\n " + str(self))





    def __str__(self):

        to_return = "ANNOTATED HIT\n"
        to_return = to_return + "Alignment start = " + str(self.five_end) + " and Alignment end: " + str(self.three_end) + "; strand:" + self.strand + "\nOperon ID:" + self.operon_id 
        to_return = to_return + super().__str__()

        return to_return
    
    def __eq__(self, other):
        if isinstance(other, AnnotatedHit):
            
            if self.genome_accession == other.genome_accession:
                return True
            else:
                return False
           
        elif isinstance(other, GenomeFeature):
            return super().__eq__(other)
            