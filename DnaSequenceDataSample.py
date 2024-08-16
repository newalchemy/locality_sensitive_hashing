# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:45 2024

@author: patri
"""
import utils;
import editdistance;
import TwoWayDict;

class DnaSequenceDataSample:
        
        
    def __init__(self, length_string, num_of_samples, shingle_size, seed=-1):
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = -1;
        self.shingle_size = shingle_size;

        
        self.sample_list = [];
        
        self.sampleID_sample_map = TwoWayDict.TwoWayDict();
        
        print('creating random samples\n');
        for i in range(0, num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sampleID_sample_map[i] = cur_sample;
            
        
        print('generated random samples\n');
        
        self.shingleThisSample()
        #print('created random samples, computing all pairs edit distance \n');
        
        #self.computeAllPairsEditDistance();
        
        #print('completed all pairs edit distance\n');
    
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
            
        
    def shingleThisSample(self):
        # Creates a dictionary of sample id - > [list of shingle ids]

        #For now, each sample can only have one shingling.  If the need arises this can be changed.
        
        #sample id to shingle id list
        self.shingling_dict = {};
        self.known_shingles = set();
        
        
        itervar = 0;
        
        self.shingleID_shingle_map = TwoWayDict.TwoWayDict();
        
        
        for i in range(0, self.num_of_samples):
            my_sample = self.getSamplebySampleID(i);
            shingles = utils.shingle_sequence(my_sample, self.shingle_size);
            
            shingle_id_set = set();
            for j in range(0, len(shingles)):
                my_shingle = shingles[j];
                
                #If we've already seen the shingle, retrieve its id for the list
                if (my_shingle in self.known_shingles):
                    shingle_id = self.getShingleIDbyShingle(my_shingle);
                    shingle_id_set.add(shingle_id);
                    
                #Otherwise, create a new id for it and add it to the known shingle list.
                else:
                    self.shingleID_shingle_map[itervar] = my_shingle;
                    shingle_id_set.add(itervar);
                
                    itervar = itervar + 1;
                    self.known_shingles.add(my_shingle);                

            shingle_id_list = list(shingle_id_set);

            update = [(i, shingle_id_list)];
            self.shingling_dict.update(update);
                            
        
        
    def getSampleIDbySample(self, sample_str):
        return self.sampleID_sample_map[sample_str];
    
    def getSamplebySampleID(self, sample_id):
        return self.sampleID_sample_map[sample_id];
    
    def getShingleIDbyShingle(self, shingle_str):
        return self.shingleID_shingle_map[shingle_str];
    
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
            
        
        
            
        