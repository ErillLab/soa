'''
Testing the function that assigns a feature to an annotated hit. 

Error: No feature found for
 ANNOTATED HIT
Alignment start = 442703 and Alignment end: 444124; strand:-
Operon ID:4941_VC1577_MCRBS_CLSTR1667
Query Accession:WP_001146752.1
Genome Feature
Vibrio_aquimaris_strain_THAF100_plasmid_pTHAF100_a_complete
NZ_CP045351.1: Coding start position: None Coding end position: None
Strand: -
Protein ID: None
Locus Tag: None

Coverage value:0.051320018780020946

'''

from Bio import Entrez

def test():

    align_five_end = 442703
    align_three_end = 444124

    #Open the record for this nucleotide record
    file_path = '/Volumes/Issac_Ext/ErillsLab/nucleotide_gbwithparts_cache/NZ_CP045351.1.xml'
    record = Entrez.read(open(file_path, 'rb'), 'xml')

    for feature in record[0]['GBSeq_feature-table']:

        if not feature['GBFeature_key'] == 'CDS':
            continue
        if not 'GBInterval_from' in feature['GBFeature_intervals'][0]:
            continue
        
        for quality in feature['GBFeature_quals']:
            if quality['GBQualifier_name'] == 'protein_id':
                protein_accession = quality['GBQualifier_value']
                if protein_accession == "WP_172971865.1":
                    print('protein found')

        coding_start = min(int(feature['GBFeature_intervals'][0]['GBInterval_from']), int(feature['GBFeature_intervals'][0]['GBInterval_to']))
        coding_end = max(int(feature['GBFeature_intervals'][0]['GBInterval_from']), int(feature['GBFeature_intervals'][0]['GBInterval_to']))
    
        #Determine which is longer: the alignment or the feature, and pulling their start and end positions. 
        long_start = align_five_end
        short_start = align_five_end

        long_end = align_three_end
        short_end = align_three_end

        if (align_three_end - align_five_end) > (coding_end - coding_start):
            print('long one is alignment')
            short_start = coding_start
            short_end = coding_end
        else:
            print('long one is feeature')
            long_start = coding_start
            long_end = coding_end


        #The feature does not include any part of the feature (or the feature does not contain any part of the alignment), skip this feature
        if short_end < long_start or short_start > long_end:
            continue

        #Calculate the coverage
        feat_coverage = -1

        #The coverage is 1 if the hit and feature line up exactly.
        #If the short one is completely included in the long one, the length of the short one divided by the length of the long one is the coverage. 
        if short_start >= long_start and short_end <= long_end:
            print('got 2')
            feat_coverage = (short_end - short_start) / (long_end - long_start)

        #If the short one "overhangs" on the five prime end, the covered region is from the begining of the long one to end of the small one. 
        elif short_start < long_start:
            print('got 3')
            feat_coverage = (short_end - long_start) / (long_end - long_start)

        #If the short one "overhangs" on the three prime end, the covered region is from the beginning of the short one to the end of the long one
        elif short_end > long_end:
            print('got 4')
            feat_coverage = (long_end - short_start) / (long_end - long_start)
        else:
            print(short_start, '\t', short_end, '\t', long_start, '\t', long_end)
        
    return feat_coverage
    
print(test())
