# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:37:23 2024

@author: patri
"""

import utils;
import editdistance;
import time;
import TwoWayDict;


test_dict = TwoWayDict();
for i in range(0, 100):
    samp = utils.generate_random_dna_sequence(100);
    test_dict[i] = samp;


dna_seq_length = 1000;

edit_dist_input = 30;

monte_carlos = 50000;
startime = time.time();

dict_log = {};
'''
print('edit distance input ', edit_dist_input, '\n');

for i in range(0, monte_carlos):
    print('mc: ', i, '\n');
    test_seq = utils.generate_random_dna_sequence(dna_seq_length);
    
    new_seq = utils.change_2N_edit_distance_in_DNA_seq(test_seq, edit_dist_input);
    
    edit_dist1 = editdistance.eval(test_seq, new_seq)
    edit_dist2 = utils.edit_distance(test_seq, new_seq);
    
    if (edit_dist1 != edit_dist2):
        print('test failed! \n');

print('test passed!');
'''

startime = time.time();

print('exceeds max diff test\n');
for i in range(0, monte_carlos):
    print('mc: ', i, '\n');
    test_seq = utils.generate_random_dna_sequence(dna_seq_length);
    new_seq = utils.generate_random_dna_sequence(dna_seq_length);
    
    edit_dist1 = editdistance.eval(test_seq, new_seq);
    
        
    
endtime = time.time();

time1 = endtime - startime;

print('elapsed time: ', time1, '\n');
print('test passed');


#print('test_seq ', test_seq, '\n');
#print('test_seq_len', len(test_seq), '\n');

#print('new_seq ', new_seq, '\n');
#print('new_seq_len', len(new_seq), '\n');

#print('edit distance ', new_edit_distance, '\n');

#print('largest shared shingle size ', shingle_size, ' \n');

#print('new string alg ran in ', elapsed, ' seconds \n');