# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:45 2024

@author: patri
"""
import utils;

class DnaSequenceDataSample:
        
        
    def __init__(self, length_string, num_of_samples, seed=-1):
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = -1;

        
        self.sample_list = [];
        
        for i in range(0, num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sample_list.append(cur_sample);
        
    
    def getLengthDNAString(self):
        return self.length_string;
    
    def getNumOfDNASamples(self):
        return self.num_of_samples;
    
    def getSeed(self):
        return self.seed;
    
    def shingleThisSample(self, shingle_size):
        # returns a dictionary {sample, shingle sequence}

        #For now, each sample can only have one shingling.  If the need arises this can be changed.
        
        self.shingle_size = shingle_size;
        self.shingling_dict = {};
        
        for i in range(0, self.num_of_samples):
            shingled = utils.shingle_sequence(self.sample_list[i], shingle_size);
            my_sample = self.sample_list[i];
            update = [(my_sample, shingled)];
            self.shingling_dict.update(update);
        
        
    def getShinglingDict(self):
        return self.shingling_dict;
    
    def getShingleSize(self):
        return self.shingle_size;
    
    def computeAllPairsEditDistance(self):
        
        return;
            
        
        
            
        