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
        self.Sample = sample;
        self.shingle_size = sample.getShingleSize();
        self.__generateLSHShingleDB();
        
        
    def getShingleSize(self):
        return self.shingle_size;
    
    def getSampleToShingleDict(self):
        return self.sample_to_shingle_dict;
    
    def getSamplesByShingle(self, shingle_str):
        shingle_id = self.sample.getShingleIDbyShingle(shingle_str);
        return self.getSamplesByShingleId(shingle_id);
    
    def getSamplesByShingleId(self, shingle_id):
        sampleIds = self.getSampleIdsByShingleId(shingle_id);
        samples = [];
        for i in range(0, len(sampleIds)):
            my_sample_id = sampleIds[i];
            sample_str = self.Sample.getSamplebySampleID(my_sample_id);
            samples.append(sample_str);
        return samples;
            
    
    def getSampleIdsByShingleId(self, shingle_id):
        sampleIds = self.shingleID_to_sampleID_dict.get(shingle_id);
        return sampleIds;

    
    def __generateLSHShingleDB(self):
        
        #Inverts the sample_to_shingle dictionry.
        #That is, instead of sample - > list of shingles, the result of this function 
        #goes shingle -> list of samples.
        
        #mapping of shingle to string
        self.shingleID_to_sampleID_dict = {-1 : -1 };
        
        
        self.num_of_sequences = self.Sample.getNumOfDNASamples();
        
        
        for i in range(0, self.num_of_sequences):
            cur_seq_id = i;
            list_of_shingle_ids = self.Sample.getShingleIdsBySampleId(cur_seq_id);
            for j in range (0, len(list_of_shingle_ids)):
                my_shingle_id = list_of_shingle_ids[j];
                
                #If the shingleid is already in the dictionary, update the dictionary.  Otherwise create a totally new entry.
                k = self.shingleID_to_sampleID_dict.get(my_shingle_id);
                if k:
                    ent = self.shingleID_to_sampleID_dict.get(my_shingle_id);
                    ent.append(cur_seq_id);
                    new_entry = [(my_shingle_id, ent)];                    
                else:
                    new_entry = [(my_shingle_id, [cur_seq_id])];
                self.shingleID_to_sampleID_dict.update(new_entry);
                
        
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
            
            
            
            
        
        
        
        
                
            
            
            
        
        
        
        
        
        