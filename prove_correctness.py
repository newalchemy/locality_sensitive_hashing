# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:58:18 2024

@author: patri
"""
import DnaSequenceDataSample;
import time;
import utils;
import pickle;
import lshShingleDB;
import collections;

length_of_dna_sequence = 1000;
num_of_dna_samples = 5000;
shingle_size = 10;


true_num_of_shingles = length_of_dna_sequence - shingle_size + 1;
seed = -1

sample = DnaSequenceDataSample.DnaSequenceDataSample(length_of_dna_sequence, num_of_dna_samples, shingle_size, seed);


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

shingleDB = lshShingleDB.lshShingleDB(sample);

print('shingleDB Generated!');

#Computationally expensive but necessary
all_shingles = list(sample.getAllShingles());
all_samples = sample.getAllSamples();


num_shingle_ids = len(all_shingles);

for i in range(0, num_shingle_ids):
    cur_shingle = all_shingles[i];
    cur_shingle_id = sample.getShingleIDbyShingle(cur_shingle);
    lsh_sample_list = shingleDB.getSamplesByShingleId(cur_shingle_id);
    
    #Verify correctness of lsh output:
    for j in range(0, len(lsh_sample_list)):
        cur_sample = lsh_sample_list[j];
        if (cur_shingle in cur_sample):
            continue;
        else:
            print('test failed in LSH return test:  Shingle: ', my_shingle, ' String: ', sample_str, ' , j= ', j , ' \n');
    
    #Now, build the list of all samples which contain the shingle via brute force and compare the lists.
    overlap_sample_list = [];
    for j in range(0, len(all_samples)):
        my_sample = all_samples[j];
        if cur_shingle in my_sample:
            overlap_sample_list.append(my_sample);
            
    passed = collections.Counter(overlap_sample_list) == collections.Counter( lsh_sample_list);
    if not passed:
        print('LSH did not return all samples');
        
print('LSH test passed');
    
        
        
    
    
    
    