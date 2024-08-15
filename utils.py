# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 09:52:23 2024



@author: patri
"""

import random;

def generate_random_dna_sequence(length, seed=-1):
    
    if (seed != -1):
        random.seed(seed);
    
    return_sequence = "";
    for i in range(0, length):
        val = random.randint(0,3);
        cur_char = int_to_dna_letter(val);
        return_sequence = return_sequence + cur_char;
    
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
   

#test

#test_seq = generate_random_dna_sequence(200);
#print(test_seq);