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

length_of_dna_sequence = 100;
num_of_dna_samples = 5000;
shingle_size = 10;


seed = 2231991

sample = DnaSequenceDataSample.DnaSequenceDataSample(length_of_dna_sequence, num_of_dna_samples, shingle_size, seed);
with open('C:/Users/patri/OneDrive/Desktop/id_mapping_backup/sample1.pkl', 'wb') as f:
    pickle.dump(sample, f);
    
print('building LSH table \n')
shingleDB = lshShingleDB.lshShingleDB(sample);
print('built LSH table \n');

with open('C:/Users/patri/OneDrive/Desktop/id_mapping_backup/lshTable_sample1.pkl', 'wb') as f:
    pickle.dump(shingleDB, f);


endtime = time.time();
