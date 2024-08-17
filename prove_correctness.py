# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:58:18 2024

@author: patri
"""
import DnaSequenceDataSampleAndLSHTable;
import time;
import utils;
import pickle;
import lshShingleDB;
import collections;

length_of_dna_sequence = 1000;
num_of_dna_samples = 50;
shingle_size = 10;


true_num_of_shingles = length_of_dna_sequence - shingle_size + 1;
seed = -1;

sample = DnaSequenceDataSampleAndLSHTable.DnaSequenceDataSampleAndLSHTable(length_of_dna_sequence, num_of_dna_samples, shingle_size, seed);


test = sample.getAllSamples();

test_seq = utils.generate_random_dna_sequence(length_of_dna_sequence);

intersect, intersect_size, shing_size = utils.find_longest_shingle_between_two_strings(test_seq, test[0]);


#Verify shingling was done correctly
for i in range(0, num_of_dna_samples):
    sample_str = sample.getSamplebySampleID(i);
    shingles = sample.getShinglesBySample(sample_str);
    num_shingles_ret = 0;

    num_of_shingles = len(shingles)
    for j in range(0, len(shingles)):
        my_shingle = shingles[j];
        if my_shingle in sample_str:
            copies = utils.occurrences_of_shingle(sample_str ,my_shingle);
            num_shingles_ret = num_shingles_ret + copies;
            continue;
        else:
            print('test failed:  Shingle: ', my_shingle, ' String: ', sample_str, ' , i= ', i, ' \n');
            raise Exception;

    if (num_shingles_ret != true_num_of_shingles):
        print('test failed! num_shingles_ret: ', num_shingles_ret, ' true_num_shingles: ', true_num_of_shingles, ' sample: ', sample_str, '\n')
        
print('test passed for data sample generation!')


#Computationally expensive but necessary

#Given any shingle in the dataset, does the LSH implementation return all document IDs associated with that shingle?
#Compare LSH answer and brute force answer to verify correctness.
all_shingles_ids = list(sample.getAllShingleIds());
all_samples = sample.getAllSamples();


num_shingle_ids = len(all_shingles_ids);

for i in range(0, num_shingle_ids):
    cur_shingle_id = all_shingles_ids[i];
    cur_shingle = sample.getShinglebyShingleID(cur_shingle_id);
    lsh_sample_list = sample.getSamplesByShingleId(cur_shingle_id);
    
    if (any(lsh_sample_list.count(x) > 1 for x in lsh_sample_list)):
        bp = 'bp';

    
    #Verify correctness of lsh output:
    for j in range(0, len(lsh_sample_list)):
        cur_sample = lsh_sample_list[j];
        try:
            if (cur_shingle in cur_sample):
                continue;
            else:
                print('test failed in LSH return test:  Shingle: ', my_shingle, ' String: ', sample_str, ' , j= ', j , ' \n');
        except TypeError:
            bp = 'bp';
    
    #Now, build the list of all samples which contain the shingle via brute force and compare the lists.
    overlap_sample_list = [];
    for j in range(0, len(all_samples)):
        my_sample = all_samples[j];
        try:
            if cur_shingle in my_sample:
                overlap_sample_list.append(my_sample);
        except TypeError:
            bp = 'bp';
    passed1 = collections.Counter(overlap_sample_list) == collections.Counter(lsh_sample_list);
    passed2 = len(lsh_sample_list) > 0;
    passed = passed1 and passed2;
    
    if not passed:
        print('LSH vs overlap test failed \n');
        print('elms in overlap: ', len(overlap_sample_list), '\n');
        print('elms in lsh list: ', len(lsh_sample_list), '\n');
        
print('LSH test passed');
    
        
        
    
    
    
    