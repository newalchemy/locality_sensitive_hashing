# -*- coding: utf-8 -*-
"""

Utils functions to be used in multiple experiments.


"""

import random;
import kshingle as ks;

def generate_random_dna_sequence(length, seed=-1):
    
    if (seed != -1):
        random.seed(seed);
    
    return_sequence_list = [];
    for i in range(0, length):
        val = random.randint(0,3);
        cur_char = int_to_dna_letter(val);
        return_sequence_list.append(cur_char);
        
    return_sequence = ''.join(return_sequence_list);
    
    return return_sequence;
    

def int_to_dna_letter(int_value):
    
    if int_value == 0:
        return 'A'
    elif int_value == 1:
        return 'T'
    elif int_value == 2:
        return 'G'
    elif int_value == 3:
        return 'C'
    else:
        raise Exception('Invalid value passed to int_to_dna_letter');

#There needs to be a more optimal way of doing this
def shingle_sequence(dna_sequence, shingle_size):
    ks_out = ks.shingleseqs_range(dna_sequence, n_min=shingle_size, n_max=shingle_size);
    out = tuple(ks_out[0]);
    return out;

#test

#test_seq = generate_random_dna_sequence(200);
#print(test_seq);