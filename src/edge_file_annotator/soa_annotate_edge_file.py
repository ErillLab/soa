'''
@author ichaudr

This script will annotate an edge file with the detailed information about the genes in the two operons. 

'''

from Bio import Entrez
from Bio.motifs import Motif, Instances
from tqdm import tqdm
from datetime import datetime
import csv
import os
import json

################ PARAMETERS ########################

#Path to csv with cluster definitions
cluster_defs_file = 'clustered_operons.csv'

#Path to saved cluster JSONs
cluster_jsons_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Saved Clusters/complete_clusters/'

#Path to save weblogos for the motifs
logo_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/src/edge_file_annotator/logos/'

#Folder with the XML genebank files for the reference genome
ref_genome_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Reference_data/el_tor_assembly/'

#Parsed genebank records for the reference genome
ref_genome_records = []

#Unannotated edges file
old_edges_file = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/Vibrio_SOA/Gephi Files/soa_network_data_1304_clusters_filtered_IC_04242021/complete_clusters_04202021_edges_avg_mins_kld.csv'

#New edges file
new_edges_file = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/src/edge_file_annotator/annotated_complete_clusters_04202021_edges_avg_mins_kld.csv'

#HTML File 
html_file_path = '/Users/ichaudr/Documents/UMBC/Lab-Erill/Isaac/issac-workspace-IE/soa/src/edge_file_annotator/annotated_complete_clusters_04202021_edges_avg_mins_kld.html'



############## HELPER FUNCTIONS ####################

def get_cluster_motifs(cluster_id, path_to_saved_clusters=cluster_jsons_path):
    '''
    Parse out the motif for a specific cluster.

    Parameters
    ----------
    cluster_id: str
        The cluster of interest.
    path_to_saved_clusters: str
        The path to the directory that is holding the saved JSONs for the saved clusters
    
    Returns
    -------
    motifs: str
        A string that is made of the concatenated consensus sequences of the motifs a part of the cluster. 
    '''

    motifs = []

    for file in os.listdir(path_to_saved_clusters):
        if file.split('.')[0] == cluster_id:
            file_reader = json.load(open(path_to_saved_clusters+file, 'r'))

            motifs_from_file = file_reader['motifs']
            
            for m in motifs_from_file:
                temp_motif = Motif(instances=Instances(m))
                #temp_motif.weblogo(fname=logo_path+cluster_id+"_"+str(motifs_from_file.index(m))+'.png')
                motifs.append(str(temp_motif.consensus))
    return ' | '.join(motifs)




def get_prot_info(prot_id, genome_records=ref_genome_records):
    '''
    This function will parse all of the meta data from the genebank record for the specified prot_id. 

    Parameters
    ----------
    prot_id: str
        The protein of interest.
    genome_records: list of XML parsed objects
        The parsed genome records for the reference species.

    Returns
    -------
    prot_info: {}
        A dictionary of all the metadata for the specified protein
    '''

    prot_info = {
        'gene':  None,
        'locus_tag': None,
        'old_locus_tag': None,
        'inference': None,
        'note': None,
        'product': None,
    }

    for record in genome_records:
        for feature in record[0]['GBSeq_feature-table']:
            for quality in feature['GBFeature_quals']:
                if quality['GBQualifier_name'] == 'protein_id':
                    if quality['GBQualifier_value'] == prot_id:
                        for quality in feature['GBFeature_quals']:
                            if quality['GBQualifier_name'] in prot_info.keys():
                                try:
                                    prot_info[quality['GBQualifier_name']] = str(quality['GBQualifier_value']).strip()
                                except:
                                    prot_info[quality['GBQualifier_name']] = None
 
    return prot_info
                            
'''
##### PARSE THE REFERENCE GENOMES ###########

for file in os.listdir(path=ref_genome_path):
    if file[-3:] == 'xml':
        ref_genome_records.append(Entrez.read(open(ref_genome_path+file, 'rb'), 'xml'))



##### GET INFO FOR EACH CLUSTER ###########

#A dictionary holding the information for each cluster.
cluster_definitions = {}

cluster_defs_reader = list(csv.reader(open(cluster_defs_file, 'r')))

#For each cluster, parse out all of the protein ids associated with it and save it the dictionary structure
for i in tqdm(range(len(cluster_defs_reader[1:])), desc='Parsing Cluster Defs'):

    line = cluster_defs_reader[i + 1]

    curr_clstr_id = line[3]

    if curr_clstr_id in cluster_definitions.keys():
        continue

    curr_prot_ids = list(line[4].replace("'", "").replace("[","").replace("]","").split(","))

    #Find the operon in this cluster with the most genes 
    for j in range((len(cluster_defs_reader[1:]))):

        if i == j:
            continue

        other_line = cluster_defs_reader[j + 1]
        other_cluster_id = other_line[3]
        other_prot_ids = list(other_line[4].replace("'", "").replace("[","").replace("]","").split(","))

        if other_cluster_id == curr_clstr_id:
            if len(other_prot_ids) > len(curr_prot_ids):
                curr_prot_ids = other_prot_ids
    
    cluster_definitions[curr_clstr_id] = {}
    for pid in curr_prot_ids:
        cluster_definitions[curr_clstr_id][pid] = get_prot_info(pid)'''

'''for c in tqdm(cluster_definitions.keys(), desc="Motif Preproccessing"):
    get_cluster_motifs(c)'''
'''
with open('cluster_defs.json', 'w') as file:
    json.dump(cluster_definitions, file)

##### RE-WRITE THE EDGES FILE ###########
old_edges_reader = list(csv.reader(open(old_edges_file, 'r')))
new_edges_writer = csv.writer(open(new_edges_file, 'a'))

for line in tqdm(old_edges_reader[1:], desc='Re-writing Edge File'):
    cluster1_info = cluster_definitions[line[0]]
    cluster1_first_prot_info = list((list(cluster1_info.values())[0]).values())

    cluster2_info = cluster_definitions[line[1]]
    cluster2_first_prot_info = list((list(cluster2_info.values())[0]).values())

    new_row = []
    for l in line:
        new_row.append(l)
    new_row.append(get_cluster_motifs(line[0]))
    new_row.append(get_cluster_motifs(line[1]))
    for l in cluster1_first_prot_info:
        new_row.append(l)
    for l in cluster2_first_prot_info:
        new_row.append(l)

    new_edges_writer.writerow(new_row)

    '''
##### WRITE THE HTML FILE WITH THE MOTIF LOGO IMAGES ###########
new_edges_file_reader = list(csv.reader(open(new_edges_file, 'r')))
html_writer = open(html_file_path, 'a')
html_content = ('<html>\n<head><!-- This file was generated by soa_annotate_edge_file.py on' +  str(datetime.now().strftime("%m/%d/%Y %H:%M:%S"))  + '-->\n<title> Annotated Edge File </title> \n</head> \n<body>\n')

#Make the table and write the header row
html_content += '\n<table class="sortable" border="1px" style="text-align: center; border: 1px solid black;width:1300; page-break-after:always; font-family: Verdana; font-size: 11; border-collapse: collapse;">\n'

html_content += '<tr>\n'

headings = 'Source,Target,Weight,Source Motifs,Target Motifs,Source Gene,Source Locus Tag,Source Old Locus Tag,Source Inference,Source Note,Source Product,Target Gene,Target Locus Tag,Target Old Locus Tag,Target Inference,Target Note,Target Product'.split(',')
for heading in headings:
    html_content += '<th>' + heading + '</th>\t'

html_content += '\n</tr>\n'


#Write the rest of the table
for row in tqdm(new_edges_file_reader, desc='Writing HTML File'):

    html_content += '<tr>\n'

    clstr1 = row[0]
    clstr2 = row[1]
    weight = row[2]

    clstr1_motif = ''
    clstr2_motif = ''

    for logo in os.listdir(path=logo_path):
        if str(clstr1+'_') in logo:
            clstr1_motif += '\t\t<a href="' + str(logo_path) + str(logo) + '">' + logo + '</a>'
            #clstr1_motif += '\t\t<img src="' + str(logo_path) + str(logo) + '" width=600 height=175 style="padding-bottom:0px; padding-right:0px;">\n'
        
        if str(clstr2+'_') in logo:
            clstr2_motif += '\t\t<a href="' + str(logo_path) + str(logo) + '">' + logo + '</a>'
            #clstr2_motif += '\t\t<img src="' + str(logo_path) + str(logo) + '" width=600 height=175 style="padding-bottom:0px; padding-right:0px;">\n'

    
    html_content += '\t<td>' + clstr1 + '</td>\n'
    html_content += '\t<td>' + clstr2 + '</td>\n'
    html_content += '\t<td>' + weight + '</td>\n'
    html_content += '\t<td>' + clstr1_motif + '</td>\n'
    html_content += '\t<td>' + clstr2_motif + '</td>\n'

    for content in row[5:]:
        html_content += '\t<td>' + str(content) + '</td>\n'


    html_content += '\n</tr>\n'

#Add the closing tags and write to file
html_content += '</table>\n</body></html>'
html_writer.write(html_content)