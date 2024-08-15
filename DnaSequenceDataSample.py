# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:45 2024

@author: patri
"""
import utils;
import editdistance;

class DnaSequenceDataSample:
        
        
    def __init__(self, length_string, num_of_samples, seed=-1):
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = -1;

        
        self.sample_list = [];
        
        for i in range(0, num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sample_list.append(cur_sample);
        
        self.computeAllPairsEditDistance();
    
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
        #Brute force
        #Compute all pairs of an edit distance.
        
        self.all_pairs_dict = {};
        list_to_update = []
        for i in range(0, self.num_of_samples):
            s1 = self.sample_list[i];
            for j in range(i + 1, self.num_of_samples):
                s2 = self.sample_list[j];
                #only because I can't sort for some reason.  will make queries take longer.  fix later.
                key = s1.join(s2);
                val = editdistance.eval(s1, s2);
                list_to_update.append((key, val));
                
        self.all_pairs_dict.update([(key, val)]);        
            
        
        
            
        