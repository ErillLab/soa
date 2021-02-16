'''
@author ichaudr

A testing module to test other modules.
'''

import soa_sim_filter
import meme_driver

sequences = [ 

    'AAAAAAAAAA',
    'CAAAAGAAAA'
]

print(soa_sim_filter.get_percent_matches(sequences[0], sequences[1]))