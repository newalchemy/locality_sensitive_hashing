# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:49:02 2024

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
length_of_dna_sequence = 1000;
num_of_dna_samples = 500000;
shingle_size = 100;
seed = 2231991

sample = DnaSequenceDataSample.DnaSequenceDataSample(length_of_dna_sequence, num_of_dna_samples, seed);
with open('C:/Users/patri/OneDrive/Desktop/id_mapping_backup/sample1_notshingled.pkl', 'wb') as f:
    pickle.dump(sample, f);

sample.shingleThisSample(shingle_size);
with open('C:/Users/patri/OneDrive/Desktop/id_mapping_backup/sample1.pkl', 'wb') as f:
    pickle.dump(sample, f);
shingleDB = lshShingleDB.lshShingleDB(sample);
with open('C:/Users/patri/OneDrive/Desktop/id_mapping_backup/lshTable_sample1.pkl', 'wb') as f:
    pickle.dump(shingleDB, f);


endtime = time.time();

elapsed = endtime - starttime;

print('The initiation took ', elapsed, ' seconds\n');

starttime = time.time();

test_seq = utils.generate_random_dna_sequence(length_of_dna_sequence);
samps = shingleDB.getLocalSamples(test_seq);

endtime = time.time();

elapsed = endtime - starttime;

print('The query took ', elapsed, ' seconds\n');
print('The query returned ', len(samps), ' samples\n')


#The initiation took  116.28058862686157  seconds

# The query took  0.09018397331237793  seconds

# The query returned  1  samples
