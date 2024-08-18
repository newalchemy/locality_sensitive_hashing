# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:37:23 2024

@author: patri
"""

import utils;
import editdistance;
import time;


dna_seq_length = 1000;

edit_dist_input = 10;

monte_carlos = 10000;
startime = time.time();

dict_log = {};
print('edit distance input ', edit_dist_input, '\n');

for i in range(0, monte_carlos):
    test_seq = utils.generate_random_dna_sequence(dna_seq_length);
    
    new_seq = utils.change_N_edit_distance_in_DNA_seq(test_seq, edit_dist_input);
    
    edit_dist1 = editdistance.eval(test_seq, new_seq);
    
    try:
        val = dict_log[edit_dist1];
    except KeyError:
        val = 0;
    
    val = val + 1;
    
    dict_log[edit_dist1] = val;
    
    #print('mc ', i, ' edit distance ', edit_dist1, '\n');
    if (len(new_seq) != len(test_seq)):
        raise Exception('test failed');
    

print(dict_log);

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