# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:30:53 2024

@author: patri
"""

#The speed/space optimization here is to use pointer tables in lieu of a direct mapping. 
#It will be done in a future update.
import utils;

class lshShingleDB:
    
    def __init__(self, sample):
        self.sample = sample;
        self.sample_to_shingle_dict = sample.getShinglingDict();
        self.shingle_size = sample.getShingleSize();
        self.__generateLSHShingleDB();
        
        
    def getShingleSize(self):
        return self.shingle_size;
    
    def getSampleToShingleDict(self):
        return self.sample_to_shingle_dict;
    
    def __generateLSHShingleDB(self):
        
        #Inverts the sample_to_shingle dictionry.
        #That is, instead of sample - > list of shingles, the result of this function 
        #goes shingle -> list of samples.
        
        #mapping of shingle to string
        self.shingle_to_sample_dict = {};
        
        self.list_of_dna_sequences = list(self.sample_to_shingle_dict.keys());
        
        self.num_of_sequences = len(self.list_of_dna_sequences);
        
        for i in range(0, self.num_of_sequences):
            cur_seq = self.list_of_dna_sequences[i];
            shingles = self.sample_to_shingle_dict.get(cur_seq);
            known_shingles = set();
            for j in range (0, len(shingles)):
                my_shingle = shingles[j];
                
                if (my_shingle in known_shingles):
                    continue;
                
                k= self.shingle_to_sample_dict.get(my_shingle);
                if k:
                    ent = self.shingle_to_sample_dict.get(my_shingle);
                    ent.append(cur_seq);
                    new_entry = [(my_shingle, ent)];                    
                else:
                    new_entry = [(my_shingle, [cur_seq])];
                self.shingle_to_sample_dict.update(new_entry);
                known_shingles.add(my_shingle);
                
        
    def getLocalSamples(self, query_dna_sequence):
        #Given a query DNA sequence, return all 'local' samples based on the LSH scheme.
        shingles = utils.shingle_sequence(query_dna_sequence, self.shingle_size);
        out_set = set();
        
        for i in range(0, len(shingles)):
            my_shingle = shingles[i];
            sample = self.shingle_to_sample_dict.get(my_shingle);
            
            out_set = set();
            if (sample != None):
                out_set.update(sample);
        
        out_samples = list(out_set);
        return out_samples;
            
            
            
            
        
        
        
        
                
            
            
            
        
        
        
        
        
        