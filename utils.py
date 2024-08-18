# -*- coding: utf-8 -*-
"""

Utils functions to be used in multiple experiments.


"""

import random;
import kshingle as ks;
from difflib import SequenceMatcher;
#from cdifflib import CSequenceMatcher;
import copy;
import math;

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
        
def dna_letter_to_int(dna_letter):
    
    if dna_letter == 'A':
        return 0;
    elif dna_letter == 'T':
        return 1
    elif dna_letter == 'G':
        return 2
    elif dna_letter == 'C':
        return 3
    else:
        raise Exception('Invalid value passed to dna_letter_to_int');
        

def change_N_edit_distance_in_DNA_seq(dna_seq, num_changes_to_make):

    out_dna_seq = copy.copy(dna_seq);
    
    seq_length = len(dna_seq);
    
    for i in range(0, num_changes_to_make):
        #change_choice_int = random.randint(0,10);
        change_choice_int = 1;
        
        inx_char_to_change = random.randint(0, seq_length - 1);
        
        #Substitution - pick a random character to substitute
        if (change_choice_int == 0):
            try:
                old_char = dna_seq[inx_char_to_change];
            except IndexError:
                print('bad index: ', inx_char_to_change, '\n');
                print('len dna seq: ', len(dna_seq));
            
            old_char_int = dna_letter_to_int(old_char);
            
            new_char_int = random.randint(0, 3);
            while (new_char_int == old_char_int):
                new_char_int = random.randint(0, 3);
            new_char = int_to_dna_letter(new_char_int);
            out_dna_seq = out_dna_seq[:inx_char_to_change] + new_char + out_dna_seq[inx_char_to_change + 1:];
            
        #Insert/Delete - insert and delete a random character.
        else:
            new_char = int_to_dna_letter(random.randint(0,3));
            inx_to_delete = random.randint(0, seq_length);
            #If the distance between the two indicies is no more than 1 then you're performing substitution, not insert/delete swap
            while (abs(inx_to_delete - inx_char_to_change) <= 1):
                inx_to_delete = random.randint(0, seq_length);
            
            if (inx_to_delete < inx_char_to_change): 
                out_dna_seq = out_dna_seq[:inx_to_delete] + out_dna_seq[inx_to_delete + 1:inx_char_to_change] + new_char + out_dna_seq[inx_char_to_change:seq_length];
            else:
                #deleting the last character is a special case
                if (inx_to_delete == seq_length):
                    out_dna_seq = out_dna_seq[:inx_char_to_change] + new_char + out_dna_seq[inx_char_to_change:seq_length - 1];
                else:
                    out_dna_seq = out_dna_seq[:inx_char_to_change] + new_char + out_dna_seq[inx_char_to_change:inx_to_delete] +  out_dna_seq[inx_to_delete + 1:seq_length];
    
    return out_dna_seq;
        
        

#There needs to be a more optimal way of doing this
def shingle_sequence(dna_sequence, shingle_size):
    ks_out = ks.shingleseqs_range(dna_sequence, n_min=shingle_size, n_max=shingle_size);
    out = tuple(ks_out[0]);
    return out;

def occurrences_of_shingle(sample, shingle):
    count = start = 0
    while True:
        start = sample.find(shingle, start) + 1
        if start > 0:
            count+=1
        else:
            return count

# REALLY inefficient but also REALLY convinent ...
def find_longest_shingle_between_two_strings(string1, string2):
    # returns intersection shingles, number of intersecting shingles, size of shingles
    intersection = 1;
    shingle_size = 1;
    prev_intersection = 0;
    prev_intersect_size = 0;
    while (True):
        shingle_size = shingle_size + 1;
        ks_out1 = set(ks.shingleseqs_range(string1, n_min=shingle_size, n_max=shingle_size)[0]);
        ks_out2 = set(ks.shingleseqs_range(string2, n_min=shingle_size, n_max=shingle_size)[0]);
        intersect = ks_out1.intersection(ks_out2);
        intersection = len(intersect);
        if (intersection == 0):
            break;
        
        prev_intersection = intersect;
        prev_intersect_size = intersection;
    
    return prev_intersection, prev_intersect_size, shingle_size;
#test

#test_seq = generate_random_dna_sequence(200);
#print(test_seq);