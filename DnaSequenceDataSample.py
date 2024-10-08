# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:45 2024

@author: patri
"""
import utils;
import editdistance;
import TwoWayDict;

class DnaSequenceDataSampleAndLSHTable:
        
        
    def __init__(self, length_string, num_of_samples, shingle_size, seed=-1):
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = -1;
        self.shingle_size = shingle_size;
        
        self.sample_list = [];
        
        self.sampleID_sample_map = TwoWayDict.TwoWayDict();
        self.generateRandomSamples();
        self.shingleThisSampleAndBuildLSHTable();
        
                #print('created random samples, computing all pairs edit distance \n');
        
        #self.computeAllPairsEditDistance();
        
        #print('completed all pairs edit distance\n');
        
    def generateRandomSamples(self):
        for i in range(0, self.num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sampleID_sample_map[i] = cur_sample;
    
    def getSampleIdSampleMap(self):
        return self.sampleID_sample_map;
    
    def getLengthDNAString(self):
        return self.length_string;
    
    def getNumOfDNASamples(self):
        return self.num_of_samples;
    
    def getSeed(self):
        return self.seed;
    
    def getAllShingles(self):
        return self.known_shingles;
    
    def getAllSamples(self):
                
        sample_list = [];
        for i in range(0, self.num_of_samples):
            sample = self.sampleID_sample_map[i];
            sample_list.append(sample);
        return sample_list;
            
        
    def shingleThisSampleAndBuildLSHTable(self):
        # Creates a dictionary of sample id - > [list of shingle ids]

        #For now, each sample can only have one shingling.  If the need arises this can be changed.
        
        #sample id to shingle id list
        self.shingling_dict = {};
        known_shingles = set();
        self.shingleID_to_sampleID_dict = {-1 : -1 };
        
        itervar = 0;
        
        self.shingleID_shingle_map = TwoWayDict.TwoWayDict();
        
        
        for i in range(0, self.num_of_samples):
            my_sample = self.getSamplebySampleID(i);
            shingles = utils.shingle_sequence(my_sample, self.shingle_size);
            
            shingle_id_set = set();
            for j in range(0, len(shingles)):
                my_shingle = shingles[j];
                
                #Create / Fetch unique shingle ID
                
                #If we've already seen the shingle, retrieve its id for the list
                if (my_shingle in known_shingles):
                    shingle_id = self.getShingleIDbyShingle(my_shingle);
                    shingle_id_set.add(shingle_id);
                    
                #Otherwise, create a new id for it and add it to the known shingle list.
                else:
                    shingle_id = itervar;
                    self.shingleID_shingle_map[shingle_id] = my_shingle;
                    shingle_id_set.add(shingle_id);
                
                    itervar = itervar + 1;
                    known_shingles.add(my_shingle);      
                    
                # Now update the shingleID to sampleId dictionary
                k = self.shingleID_to_sampleID_dict.get(shingle_id);
                if k:
                    ent = self.shingleID_to_sampleID_dict.get(shingle_id);
                    ent.append(my_sample);
                    new_entry = [(shingle_id, ent)];                    
                else:
                    new_entry = [(shingle_id, [i])];
                self.shingleID_to_sampleID_dict.update(new_entry);


            shingle_id_list = list(shingle_id_set);

            update = [(i, shingle_id_list)];
            self.shingling_dict.update(update);
                            
        
        
    def getSampleIDbySample(self, sample_str):
        return self.sampleID_sample_map[sample_str];
    
    def getSamplebySampleID(self, sample_id):
        return self.sampleID_sample_map[sample_id];
    
    def getShingleIDbyShingle(self, shingle_str):
      try:
        return self.shingleID_shingle_map[shingle_str];
      except KeyError:
        return -1;
    
    def getShinglebyShingleID(self, shingle_id):
        return self.shingleID_shingle_map[shingle_id];
    
    def getShingleIdsBySampleId(self, sample_id):
        return self.shingling_dict[sample_id];

    def getShinglesBySampleId(self, sample_id):
        shingle_id_list = self.shingling_dict[sample_id];
        shingles_list = [];
        for j in range(0, len(shingle_id_list)):
            my_shingle_id = shingle_id_list[j];
            my_shingle = self.shingleID_shingle_map[my_shingle_id];
            shingles_list.append(my_shingle);
        return shingles_list;
    
    def getShinglesBySample(self, sample_str):
        sample_id = self.getSampleIDbySample(sample_str);
        return self.getShinglesBySampleId(sample_id);

        
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
            
        
        
            
        