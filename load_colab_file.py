# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:39:50 2024

@author: patri
"""

import pickle
import lshShingleDB
from DnaSequenceDataSample import DnaSequenceDataSample;
from TwoWayDict import TwoWayDict;
from lshShingleDB import lshShingleDB;
import utils;

path_sample = 'C:/Users/patri/OneDrive/Desktop/LSH_Data/sample_3/lshTable_sample3.pkl';
path_lsh = 'C:/Users/patri/OneDrive/Desktop/LSH_Data/sample_3/lshTable_sample3.pkl';

#file = open(path_sample, 'rb');

#Sample_Obj = pickle.load(file)

objects = []
with (open(path_sample, "rb")) as openfile:
    while True:
        try:
            objects.append(pickle.load(openfile))
        except EOFError:
            break

#file.close();

file = open(path_lsh, 'rb');

Lsh_Obj = pickle.load(file)

file.close();

#Lsh_Obj.setSampleTable(Sample_Obj);

#length_of_dna_sequence = Sample_Obj.getLengthDNAString();

#test_seq = utils.generate_random_dna_sequence(length_of_dna_sequence);
#samps = Lsh_Obj.getLocalSamples(test_seq);

bp = 'bp';
