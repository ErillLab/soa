'''
@author ichaudr

A systematic analysis of operon regulation. A set of experimentally supported and computationally predicted set of operons are passed in and used to determine sets of homologous operons (via BLAST).

For each set of homologous operon group, the promoter region of each leading operon gene will be used in motif discovery.

The pairwise distance between all motifs are calculated and used to plot operons in a network where the distance between two operons is proportional to the pairwise distance between their motifs. 

'''

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline, NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Seq, Entrez
from soa_features import AnnotatedHit
from soa_genome_fragment import GenomeFragment
from soa_sim_filter import sim_filter
from soa_operon_cluster import OperonCluster
from soa_motif_analysis import *
from tqdm import tqdm
import time
import re
import csv
import os
import os.path
import pickle
import json
import sys


#Entrez Parameters
Entrez.email = 'ichaudr1@umbc.edu'
Entrez.api_key = 'cd2f5bf7a67d086647ec33da2c985e018d08'
REQUEST_LIMIT = 5
SLEEP_TIME = 0.5

#Reference Parameters
reference_assembly_accession = 'GCF_000006745.1'
reference_genome_accessions = ['NC_002505.1', 'NC_002506.1']
reference_genome_name = 'NamePlaceHolder'
reference_operons_file = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/reference_operons/test.csv'

#BLAST Parameters
local_blast_root = '../local_blast_bin/'
blast_results_filename = local_blast_root + '{query_accession}_{db_id}_{e_value}_{coverage}.p'
minimum_coverage = 0.75
blastdb_path = '/Volumes/Issac_Ext/ErillsLab/SOA/vibrio_nucl/vibrio_database_nucl'
blast_num_threads = 5

#Reverse BLAST parameters
rev_blast_root = '../rev_blast_bin/{accession}/'
rev_blast_all_protein_file = rev_blast_root + 'seqs.fasta'
rev_blast_db = rev_blast_root + '{accession}_db'
rev_blast_input_file = rev_blast_root + 'input.fasta'
rev_blast_output_file = rev_blast_root + 'output'
rev_blast_num_threads = 5

#Operon Assembly and Promoter Parameters
operon_min_width = 10
operon_max_width = 300
max_seq_sim = .85

#MEME and motif parameters
meme_exec_path = ''
motif_e_val_threshold = 10E-3

#Output Parameters
cache_dir = '/Volumes/Issac_Ext/ErillsLab/nucleotide_gbwithparts_cache/'
promoter_out_dir = '../output/promoter_seqs/{id}_promoters.fasta'
meme_out_dir = '../meme_bin/{id}_meme_out/'


def local_blast_search(input_record, db_path, operon_id, db_id='Vibrio', e_cutoff=10E-10, min_cover=None):
    '''
    Completes a local blast search instead of conducting remotely. Includes a reverse BLAST filter to remove any hits that are not homologs. 

    Parameters
    ----------
    input_record: str
        The accession for the query protein record
    db_path: str
        The file path to the local BLAST database
    operon_id: str
        The ID for the reference operon the gene belongs to
    e_cutoff: float
        The maximum e-value to cutoff the blast hits.
    min_cover: float
        The minimum coverage for each accepted hit
    
    Returns
    -------
    annotated_hits: list[AnnotateHit.object]
        A list of AnnotatedHit.objects that hold metadata for each of the BLAST hits.
    '''

    #Formatted BLAST results file
    f_blast_results_file = blast_results_filename.format(query_accession=input_record, db_id=db_id, e_value=str(e_cutoff), coverage=str(min_cover))

    #Check if the BLAST results have been saved from a previous run. If so, return the saved results.
    if os.path.isfile(f_blast_results_file):
        print('Loading saved results: ' + str(f_blast_results_file))
        records = pickle.load(open(f_blast_results_file, 'rb'))
        print('--- Returning ' + str(len(records)) + ' saved hits for ' + str(input_record))
        return records

    print('Checking if reference db is made...')
    make_reference_blastdb()

    print('Local BLAST: ' + str(input_record))

    #Check if the FASTA record for the input record is saved locally
    query_file = local_blast_root + input_record + '_seq.fasta'
    if not os.path.isfile(query_file):

        #Get the FASTA record for the input record
        fasta_record = None

        for i in range(REQUEST_LIMIT):

            try:

                handle = Entrez.efetch(db='protein', id=input_record, retmode='fasta', rettype='xml')        
                fasta_record = handle.read()
                time.sleep(SLEEP_TIME)

                break

            except:

                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

                if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        #Check if the fasta record was pulled successfully
        if fasta_record == None:
            print('\t\tFasta record could not be downloaded for ' + str(input_record))
            return None
        
        #Write the FASTA record to a temporary input file for the local blast search
        with open(query_file, 'w') as file:
            file.write(fasta_record)

        print('Downloaded query sequence: ' + str(input_record))

    #Get the query length that will be used to calculate the coverage later
    record_fasta = SeqIO.read(open(query_file,'r'),'fasta')
    query_length = len(record_fasta.seq)

    #Conduct the local BLAST search
    blast_outfile = local_blast_root + 'out.xml'
    blast_command = NcbitblastnCommandline(query=query_file, db=db_path, evalue=e_cutoff, outfmt=5, out=blast_outfile, num_threads=blast_num_threads)
    blast_command()

    #Parse the BLAST results
    blast_records = list(NCBIXML.parse(open(blast_outfile, 'r')))

    print('Parsing through local BLAST results ' + str(input_record) + '...')

    #List of annotated hits to return
    return_hits = []

    print("\t|~> Extracting hits from BLAST results...")
    for record in blast_records[0].alignments:
        current_hit_def = re.sub('[^A-Za-z0-9]+', '_', '_'.join(record.hit_def.split(' ')[1:-1]))
        curr_hit_rec = record.hit_def.split(' ')[0]
        print("\t\t|~> Analyzing hit " + str(curr_hit_rec) + ' in ' + str(current_hit_def))
        
        #Iterate through the hits
        for hit in record.hsps:

            #Initiates a AnnotatedHit object if set by the parameters.
            a_hit = AnnotatedHit(query_accession=input_record, operon_id=operon_id, hit_accession=curr_hit_rec, genome_fragment_name=current_hit_def, align_start=hit.sbjct_start, alignment_seq=hit.sbjct, 
                align_end=hit.sbjct_end, strand=hit.frame[1], percent_identity=(hit.identities/hit.align_length), req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME)
            
            if min_cover == None:
                continue

            #Calculate the coverage for the current hit                  
            cov = (hit.query_end - hit.query_start + 1) / (query_length)
            print('\t\t\tCoverage value: ' + str(cov))
            
            if(cov >= min_cover):
                
                return_hits.append(a_hit)
                    
            #Prints error if the minimum coverage is not met    
            else:
                print("\t\t|~> Hit did not meet coverage requirement: " + str(curr_hit_rec))
                print('\t\t\tCoverage value: ' + str(cov))
    
    #Reverse BLAST filter
    purged_hits = []
    for ahit in tqdm(return_hits, desc=('Reverse BLAST ' + input_record)):
        if check_reverse_blast(input_record, ahit):
            purged_hits.append(ahit)
    
    print('\t|~> Reverse BLAST filtered out ' + str(len(return_hits) - len(purged_hits)) + ' hits.')

    return_hits = purged_hits

    #Save the BLAST results to a file
    pickle.dump(return_hits, open(f_blast_results_file, 'wb'))

    print("\t|~> Returning " + str(len(return_hits)) + " unique hits")
    return return_hits


def make_reference_blastdb():
    '''
    Makes the local BLAST database out of the nucleotide records associated with the reference organism. This will be used for the reverse BLAST.
    '''
    '''
    global reference_assembly_accession
    global rev_blast_root
    global rev_blast_db
    global rev_blast_all_protein_file'''

    #Check if a local database has already been made for this reference. Exit function if it already exits
    if os.path.exists(rev_blast_root.format(accession=reference_assembly_accession)):
        return 
    
    #Make the directory for the BLAST db
    os.mkdir(rev_blast_root.format(accession=reference_assembly_accession))

    #Gets all genomes associated with the reference genome assembly by searching the nucleotide DB with the genome assembly accession
    for i in range(REQUEST_LIMIT):

        try:
            handle = Entrez.esearch(db='nuccore', term=reference_assembly_accession, retmode='xml', retmax='5000')
            search_records = Entrez.read(handle)
            time.sleep(SLEEP_TIME)
        except:
            print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

            if i == ( REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")
    
    #Download all of the protein sequences from each of the nucelotide records associated with this assembly and write them to a single fasta file to be used to compile the database. 
    #Stores all the sequence objected for all the proteins in the species
    sequences_to_write = []

    #Fetches the records
    for i in range(REQUEST_LIMIT):

        try:
            handle = Entrez.efetch(db='nuccore', id=search_records['IdList'], start='begin', stop='end', rettype='fasta_cds_aa')
            records = list(SeqIO.parse(handle, 'fasta'))
            time.sleep(SLEEP_TIME)
            break
        except:
            print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

            if i == ( REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")

    sequences_to_write.extend(records)

    #Write the sequences to the output FASTA file
    if len(sequences_to_write) > 0:
        print('\tWriting all CDS features for ' + str(reference_assembly_accession))
        SeqIO.write(sequences_to_write, rev_blast_all_protein_file.format(accession=reference_assembly_accession), 'fasta')
    else:
        print('\tNo sequences found for ' + str(reference_assembly_accession))

    #Make the local BLAST db
    mkblastdb_command = NcbimakeblastdbCommandline(dbtype='prot', input_file=rev_blast_all_protein_file.format(accession=reference_assembly_accession), out=rev_blast_db.format(accession=reference_assembly_accession))
    mkblastdb_command()

    print('\tMade local db for ' + str(reference_assembly_accession))


def check_reverse_blast(query_accession, annotated_hit):
    '''
    Performs a local BLAST search of the hit against the query to see if the same feature is returned. 

    Parameters
    ----------
    query_accession: string
        The accession for the query gene that the hit resulted from.
    annotated_hit: AnnotatedHit object
        The hit that is being tested. 
    '''
    '''
    global rev_blast_db
    global reference_assembly_accession
    global rev_blast_input_file
    global rev_blast_output_file'''

    #Input and output files
    temp_in = rev_blast_input_file.format(accession=reference_assembly_accession)
    temp_out = rev_blast_output_file.format(accession=reference_assembly_accession)

    #Clean the alignment sequence of the '-' from the gaps. This avoids a warning when running the local BLAST
    purged_alignment_seq = ''

    for c in annotated_hit.alignment_seq:
        if not c == '-':
            purged_alignment_seq = purged_alignment_seq + c

    #Write the sequence of the hit sequence to the temp_in file
    SeqIO.write(SeqRecord(Seq.Seq(purged_alignment_seq), id='temp'), temp_in, 'fasta')

    #Path to the database
    db_path = rev_blast_db.format(accession=reference_assembly_accession)

    #Complete BLAST search
    revblast_cmd = NcbiblastpCommandline(query=temp_in, db=db_path, out=temp_out, outfmt=5, num_threads=rev_blast_num_threads)
    revblast_cmd()

    #Open the output file and determine if there were any hits to the query
    hit_returned_query = False
    with open(temp_out) as file:

        #Parse the output file from the BLAST search
        result_handle = list(NCBIXML.parse(file))

        if len(result_handle[0].alignments) > 0:
            record  = result_handle[0].alignments[0]
            if (query_accession in record.hit_def):
                hit_returned_query = True
                
    return hit_returned_query

def get_reference_intergenic_distace(genes):
    '''
    Determines the average intergeneic distance between a set of genes in an operon. 

    Parameters
    ----------
    genes: [str]
        The list of protein accessions that represent the products of the genes of interest.
    
    Returns
    -------
    max_intergenic_distance: int
        The maximum intergenic distance.
    '''

    #A list of GenomeFragment objects for each of the genome accessions associated with the reference
    reference_genome_frags = []

    for accession in reference_genome_accessions:
        frag = GenomeFragment(name=reference_genome_name, genome_fragment_accession=accession, req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME, cache_directory=cache_dir)
        frag.fetch_features()
        reference_genome_frags.append(frag)

    if len(reference_genome_frags) == 0:
        print('No reference GenomeFragments made.')
        return -1
    
    #Determine which reference GenomeFragment contains the reference operon of interest
    target_frags = []
    
    for frag in reference_genome_frags:
        
        #Keep track if all the reference genes are present in this fragment
        all_genes_present = True

        #Determine if all the reference genes are present in this fragment
        for gene in genes:
            
            gene_found = False
            
            for feat in frag.all_features:
                if feat.protein_accession == gene:
                    gene_found = True

            if not gene_found:
                all_genes_present = False

        #If all the genes were found in the fragment, it is added to the list of target fragments
        if all_genes_present:
            target_frags.append(frag)
        
    #Holds the maximum intergenic distance
    max_intergenic_distance = -1

    #Determine each of the intergenic distances in each of the target fragments
    for frag in target_frags:

        #Pull the reference features in this fragment
        ref_feats = []

        for gene in genes:
            for feat in frag.all_features:
                if gene == feat.protein_accession:
                    ref_feats.append(feat)
        
        if not len(ref_feats) == len(genes):
            print("ABCDEF")

        for i in range(len(ref_feats) - 1):
            temp_distance = ref_feats[i].get_intergenic_distance(ref_feats[i+1])
            if temp_distance > max_intergenic_distance:
                max_intergenic_distance = temp_distance

    return max_intergenic_distance

def load_input_json(path):
    '''
    Loads the parameters from the input JSON file.

    Parameters
    ---------
    path: str
        Path to the input JSON
    
    Returns
    -------
    None
    '''

    file_reader = json.load(open(path))

    global REQUEST_LIMIT
    global SLEEP_TIME
    REQUEST_LIMIT = file_reader['entrez'][0]['request_limit']
    SLEEP_TIME = file_reader['entrez'][0]['sleep_time']
    Entrez.email = file_reader['entrez'][0]['email']
    Entrez.api = file_reader['entrez'][0]['api_key']

    global minimum_coverage
    global blastdb_path
    global blast_num_threads
    global rev_blast_num_threads
    minimum_coverage = file_reader['blast'][0]['minimum_coverage']
    blastdb_path = file_reader['blast'][0]['blastdb_path']
    blast_num_threads = file_reader['blast'][0]['blast_num_threads']
    rev_blast_num_threads = file_reader['blast'][0]['rev_blast_num_threads']

    global reference_assembly_accession
    global reference_genome_accessions
    global reference_genome_name
    global reference_operons_file
    reference_assembly_accession = file_reader['reference_species'][0]['reference_genome_assembly_accession']
    reference_genome_accessions = file_reader['reference_species'][0]['reference_genome_accessions']
    reference_genome_name = file_reader['reference_species'][0]['reference_genome_name']
    reference_operons_file = file_reader['reference_species'][0]['reference_operons_file']

    global operon_min_width
    global operon_max_width
    global max_seq_sim
    operon_min_width = file_reader['operon_prom_assemb'][0]['operon_min_width']
    operon_max_width = file_reader['operon_prom_assemb'][0]['operon_max_width']
    max_seq_sim = file_reader['operon_prom_assemb'][0]['max_seq_sim']

    global meme_exec_path
    global motif_e_val_threshold
    meme_exec_path = file_reader['meme_motif'][0]['meme_exec_path']
    motif_e_val_threshold = file_reader['meme_motif'][0]['motif_e_val_threshold']

    global cache_dir
    cache_dir = file_reader['output'][0]['cache_dir']


def soa():
    '''
    The main function in sys_analysis.py
    
    Pipeline:
    1. Read in operons from the input CSV file.
    2. For each gene in each operon, perform a BLAST search.
    3. Confirm BLAST results via reverse BLAST. 
    4. Group all BLAST hits based off of nucleotide accession (GenomeFragment objects).
    5. For each GenomeFragment, fetch all features and assemble operons.
    '''

    #A list of operon for all of the reference operons in the following format:
    # {'operon_id': id, 'gene_accessions': [accessions]}
    all_operon_info = []

    #A list of all the BLAST search results
    blast_hits = []
    
    '''
    #The list of promoters found for each of the homologous sets for each reference operon. 
    #Format {'operon_id':[list of promoters]}
    indexed_promoters = {}'''

    #The list of operon clusters that contains all of the homologous operons in each cluster
    operon_clusters = []

    #Complete a BLAST search for each of the genes in each of the reference operons 
    with open(reference_operons_file, 'r') as file:
        
        csv_reader = csv.reader(file)

        for operon_data in csv_reader:

            #Parse out operon id
            operon_id = operon_data[0]

            #Parse out gene ids for this operon
            gene_ids = operon_data[2].replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')

            operon_info = {
                'operon_id': operon_id,
                'gene_accessions': gene_ids,
                'intergenic_distance': 100
            }

            all_operon_info.append(operon_info)
            
            for gene_id in gene_ids:
                blast_hits.extend(local_blast_search(input_record=gene_id, db_path=blastdb_path, operon_id=operon_id, min_cover=minimum_coverage))

    print('Total number of hits returned: ' + str(len(blast_hits)))

    #Calculate all intergenic distances
    for op_info in tqdm(all_operon_info, desc='Calculating Intergenic Distance'):
        op_info['intergenic_distance'] = get_reference_intergenic_distace(genes=op_info['gene_accessions'])

    #A list of all the GenomeFragments
    genome_fragments = []

    #Sort all of the BLAST hits into GenomeFragments
    for hit in tqdm(blast_hits, desc='Sorting hits into Genome Fragments'):

        #Keeps track of whether or not the hit has been added
        hit_added = False

        #Check if there is a GenomeFragment with the same nucleotide accession as the hit

        if len(genome_fragments) == 0:
            new_genome_fragment = GenomeFragment(name=hit.genome_fragment_name, genome_fragment_accession=hit.genome_accession, req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME, cache_directory=cache_dir)
            new_genome_fragment.add_hit(hit)
            genome_fragments.append(new_genome_fragment)
            hit_added = True

        for frag in genome_fragments:

            if hit_added:
                continue
            
            if frag.genome_accession == hit.genome_accession:
                frag.add_hit(hit)
                hit_added = True
            
        if not hit_added:
            new_genome_fragment = GenomeFragment(name=hit.genome_fragment_name, genome_fragment_accession=hit.genome_accession, req_limit=REQUEST_LIMIT, sleep_time=SLEEP_TIME, cache_directory=cache_dir)
            new_genome_fragment.add_hit(hit)
            genome_fragments.append(new_genome_fragment)
            hit_added = True
    
    print('Total number of fragments: ' + str(len(genome_fragments)))

    #Process each of the fragments
    for frag in tqdm(genome_fragments, desc='Processing Fragments'):

        tqdm.write('Current Fragment: ' + str(frag.name) + '|' + str(frag.genome_accession))
        
        #Fetch genome features
        tqdm.write('----- Fetching features')
        frag.fetch_features()
        frag.fetch_hit_features()

        #Sort the hits in the fragment by operon_id
        frag.sort_hits()

        #Iterate through all of the operon_ids and assemble each operon
        tqdm.write('----- Assembling operons')
        for op_info in all_operon_info:
            num_added = frag.assemble_operons(operon_id=op_info['operon_id'], feature_limit=3, intergenic_limit=op_info['intergenic_distance'])
            tqdm.write('---------- Operon ID: ' + str(op_info['operon_id']) + '--> ' + str(num_added) + ' operons.')
        
        #Fetch the promoters for all the operons
        tqdm.write('----- Fetching operon promoters')
        frag.get_promoters()
        for operon in frag.operons:
            
            if operon.promoter == None:
                continue
            
            
            if len(operon.promoter) < operon_min_width or len(operon.promoter) > operon_max_width:
                continue

            op_added_to_cluster = False

            for cluster in operon_clusters:
                if cluster.cluster_id == operon.cluster_id:
                    cluster.operons.append(operon)
                    op_added_to_cluster = True
            
            if not op_added_to_cluster:
                new_cluster = OperonCluster(cluster_id=operon.cluster_id)
                new_cluster.operons.append(operon)
                operon_clusters.append(new_cluster)
                op_added_to_cluster = True
                
            '''
            if not operon.operon_id in indexed_promoters.keys():
                indexed_promoters.update({operon.operon_id:[]})

            indexed_promoters[operon.operon_id].append(str(operon.promoter))'''

        tqdm.write('__________'*10)
        frag.clean()


    #Filter promoters for each operon cluster
    for cluster in tqdm(operon_clusters, desc='Filtering promoters'):
        tqdm.write('Filtering ' + cluster.cluster_id)
        cluster.filter_promoters(threhold_percent_id=max_seq_sim)
        cluster.write_promoters(output_file=promoter_out_dir.format(id=cluster.cluster_id))
    
    #Preform MEME analysis on all of the clusters
    for cluster in tqdm(operon_clusters, desc='MEME analysis'):
        tqdm.write('MEME: ' + str(cluster.cluster_id))

        #Setup parameters
        in_file = promoter_out_dir.format(id=cluster.cluster_id)
        out_dir = meme_out_dir.format(id=cluster.cluster_id)

        run_meme(input_file=in_file, output_dir=out_dir)
        cluster.motifs = get_motifs(meme_data_dir=out_dir, e_val_threshold=motif_e_val_threshold)

    for c in operon_clusters:
        cluster_file_name = '../output/complete_clusters/{cluster_id}.p'
        pickle.dump(c, open(cluster_file_name.format(cluster_id=c.cluster_id)))
        print(c)
        print('~'*15)
        

    '''
    for k in indexed_promoters:
        print(k + ": " + str(len(indexed_promoters[k])))

    #Filter each list of promoters for each operon
    for operon_id in tqdm(indexed_promoters.keys(), desc='Filtering promoters'):
        filtered_list = sim_filter(seqeunces=indexed_promoters[operon_id], threhold_percent_id=0.85)
        indexed_promoters[operon_id] = filtered_list
    
    for k in indexed_promoters:
        print(k + ": " + str(len(indexed_promoters[k])))

    #Write a FASTA file for each operon containing all of the promoter sequences for that operons.
    print('Writing promoter sequnces to FASTA...')
    for operon_id in indexed_promoters:
        print('----- ' + str(operon_id))
        all_recs = []
        i = 0
        for prom in indexed_promoters[operon_id]:
            all_recs.append(SeqRecord(seq=Seq.Seq(prom), id=(operon_id + '_' + str(i)), description='|'))
            i += 1
        SeqIO.write(all_recs, promoter_out_dir.format(operon_id=operon_id), 'fasta')
        all_recs = [] '''


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print('Incorrect input. Usage: sys_analysis.py [inputfile_name].json')
        quit

    full_file_path = '../input/' + sys.argv[1]

    print('Loading input file "', full_file_path, '" ...')
    load_input_json(full_file_path)
    soa()
