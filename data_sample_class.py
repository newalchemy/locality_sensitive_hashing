# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:45 2024

@author: patri
"""
import utils;

class DnaSequenceDataSample:
        
    
    def __init__ (self, length_string, num_of_samples, seed=-1):
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = -1;
        
        self.sample_list = [];
        
        for i in range(0, length_string):
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
        
        # Note:  Future improvement to this may be to store all returns each time this method is called, 
        # and provide a separate getter to return all shinglings.
        
        ret = {};
        
        for i in range(0, self.num_of_samples):
            shingled = utils.shingle_sequence(self.sample_list[i], shingle_size);
            ret.update(self.sample_list[i], shingled);
        
        return ret;
    
            
        
        
            
        