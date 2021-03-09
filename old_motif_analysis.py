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
    
    Returns
    -------
    max_ic: float
        The informtion content for the alignment
    '''

    #The sequences of motif_a and motif_b that are in the alignment given by the offset 
    a_seqs = []
    b_seqs = []

    if offset < 0:
        offset = -offset
        alignment_length = min(len(motif_b)-offset, len(motif_a))
        a_seqs = [seq[:alignment_length] for seq in motif_a.instances]
        b_seqs = [seq[offset:alignment_length+offset] for seq in motif_b.instances]
    else:
        alignment_length = min(len(motif_a)-offset, len(motif_b))
        a_seqs = [seq[offset:alignment_length+offset] for seq in motif_a.instances]
        b_seqs = [seq[:alignment_length] for seq in motif_a.instances]
    
    #Create a temp Motif object with the combined sequences from motif_a and motif_b. 
    temp_motif = Motif(instances=Instances(a_seqs+b_seqs))
    temp_motif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)
    #return temp_motif.pssm.mean()

    return abs(sum([val for vals in temp_motif.pssm.values() for val in vals if type(val) is float and abs(val) >= 0]))

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
    
    #Dictionary in the format: {offset_value:ic},  where ic is the information content at that offset value.
    offset_ic = {}

    for offset in range(-len(motif_b) + 1, len(motif_a)):
        
        #Calculate the information content for the alignment with the current offset.
        offset_ic[offset] = get_ic(motif_a, motif_b, offset)

    for o in offset_ic:
        print(o, '\t', offset_ic[o])

    #Return all offsets that have the maximum in their alignment
    return [offset for offset in offset_ic.keys() if offset_ic[offset] == max(offset_ic.values())]



def calculate_motif_distance(motif_a, motif_b):
    '''
    Returns the distance between two motifs by:
        1. Finding the optimal alignment that maximizes the information content. 
        2. Computes the Euclidian distance between the aligned columns. 
    
    Parameters
    ----------
    motif_one, motif_two: Motif objects

    Returns
    -------
    motif_distance: float
        The computed distance between the two motifs.
    '''

    #Determine the optimal alignments for the two motifs
    alignment_offsets = get_alignment_offset(motif_a, motif_b)

    if len(alignment_offsets) == 0:
        print('Error - no alignemnt was found.')
        return -1
    
    #Holds all calculated distances
    distances = []

    for offset in alignment_offsets:
        print("Current offset: ", offset)
        for pos in range(min(len(motif_a), len(motif_b) - offset)):
            cola = dict((let, motif_a.pwm[let][pos]) for let in "ACTG")
            colb = dict((let, motif_b.pwm[let][pos + offset]) for let in "ACTG")

            distances.append(math.sqrt(sum((cola[let]-colb[let])**2 for let in "ACTG")))

    return sum(distances) / len(distances)


def calculate_motif_distance_padded(motif_a, motif_b):
    '''
    Returns the distance between two motifs by:
        1. Finding the optimal alignment that maximizes the information content. 
        2. Computes the Euclidian distance between the aligned columns. 
    
    Parameters
    ----------
    motif_one, motif_two: Motif objects

    Returns
    -------
    motif_distance: float
        The computed distance between the two motifs.
    '''

    #Determine the optimal alignments for the two motifs
    alignment_offsets = get_alignment_offset(motif_a, motif_b)

    if len(alignment_offsets) == 0:
        print('Error - no alignemnt was found.')
        return -1
    
    #Make padded versions of the motifs that take into account the offset.
    #For example: TGGAGGGCTT------------
    #             ---------TTTAGGATCTCTC
    #Would be padded as follows before the distance is calculated:
    #   TGGAGGGCTTNNNNNNNNNNNN 
    #   NNNNNNNNNTTTAGGATCTCTC
    #
    #Make a new pwm (PositionalWeightMatrix) for each matrix to account for equal probabilities 
    





    #Holds all calculated distances
    distances = []

    for offset in alignment_offsets:
        print("Current offset: ", offset)
        for pos in range(min(len(motif_a), len(motif_b) - offset)):
            cola = dict((let, motif_a.pwm[let][pos]) for let in "ACTG")
            colb = dict((let, motif_b.pwm[let][pos + offset]) for let in "ACTG")

            distances.append(math.sqrt(sum((cola[let]-colb[let])**2 for let in "ACTG")))

    return sum(distances) / len(distances)

########################################################################################################
