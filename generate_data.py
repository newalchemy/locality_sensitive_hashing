# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 16:50:56 2024

@author: patri
"""

import lshShingleDB;
import DnaSequenceDataSample;
import time;
import utils;
import pickle;

starttime = time.time();

#assume a max edit distance of 49
#TheN S = N/(R+1) = 120
#length_of_dna_sequence = 1000;
#num_of_dna_samples = 500000;
#shingle_size = 100;

#length_of_dna_sequence = 1000;
#num_of_dna_samples = 500000;
#shingle_size = 91;

length_of_dna_sequence = 1000;
num_of_dna_samples = 5000;
shingle_size = 10;


seed = 2231991

storage_path = 'C:/Users/patri/OneDrive/Desktop/LSH_Data/';
sample_name = 'test1'

sample = DnaSequenceDataSample.DnaSequenceDataSample(length_of_dna_sequence, num_of_dna_samples, shingle_size, seed);

all_samples = sample.getAllSamples();
all_lengths = [];

test_seq = utils.generate_random_dna_sequence(num_of_dna_samples);

for i in range(0, len(all_samples)):
    my_samp = all_samples[i];
    shing = utils.find_longest_shingle_between_two_strings(test_seq, my_samp);
    lshing = len(shing);
    all_lengths.append(lshing);
    
print(all_lengths);
    



strname = storage_path + sample_name + "_dna_sample_table.pkl";

with open(strname, 'wb') as f:
    pickle.dump(strname, f);

with open(strname, 'rb') as f:
    pickle.load(strname, f)

with open('C:/Users/patri/OneDrive/Desktop/LSH_Data/sample1.pkl', 'wb') as f:
    pickle.dump(sample, f);

print('building LSH table \n')
shingleDB = lshShingleDB.lshShingleDB(sample);
print('built LSH table \n');

with open('C:/Users/patri/OneDrive/Desktop/LSH_Data/lshTable_sample1.pkl', 'wb') as f:
    pickle.dump(shingleDB, f);


endtime = time.time();

elapsed = endtime - startime;

print('the program finished in ', elapsed, 'seconds ');