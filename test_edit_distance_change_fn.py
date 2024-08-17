# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:37:23 2024

@author: patri
"""

import utils;
import editdistance;
from fixed_edit_distance import randomWithFixedDistance;
import time;


dna_seq_length = 1000;

monte_carlos = 100000;
startime = time.time();


for i in range(0, monte_carlos):
    test_seq = utils.generate_random_dna_sequence(dna_seq_length);


    new_seq = utils.change_N_edit_distance_in_DNA_seq(test_seq, 7);

    if (len(new_seq) != len(test_seq)):
        raise Exception('test failed');
    

endtime = time.time();
elapsed = endtime - startime;

print('test passed');


#print('test_seq ', test_seq, '\n');
#print('test_seq_len', len(test_seq), '\n');

#print('new_seq ', new_seq, '\n');
#print('new_seq_len', len(new_seq), '\n');

#print('edit distance ', new_edit_distance, '\n');

#print('largest shared shingle size ', shingle_size, ' \n');

#print('new string alg ran in ', elapsed, ' seconds \n');