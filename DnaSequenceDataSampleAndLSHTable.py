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
        
        
        self.sampleID_sample_map = TwoWayDict.TwoWayDict();
        self.shingleID_shingle_map = TwoWayDict.TwoWayDict();
        self.sampleID_to_shingleID_dict = {};
        #this is the LSH table.  Needs an initial entry to avoid key error.
        self.shingleID_to_sampleID_dict = {-1 : -1 };


        self.__generateRandomSamples();
        self.__shingleThisSampleAndBuildLSHTable();
        
                #print('created random samples, computing all pairs edit distance \n');
        
        #self.computeAllPairsEditDistance();
        
        #print('completed all pairs edit distance\n');
        
    def getLocalSamples(self, query_dna_sequence):
        #Given a query DNA sequence, return all 'local' samples based on the LSH scheme.
        shingles = utils.shingle_sequence(query_dna_sequence, self.shingle_size);
        out_set = set();
        
        for i in range(0, len(shingles)):
            my_shingle = shingles[i];
            sample_list = self.getSamplesByShingle.get(my_shingle);
            
            if (sample_list != None):
                out_set.update(query_dna_sequence);
        
        out_samples = list(out_set);
        return out_samples;


    def __generateRandomSamples(self):
        for i in range(0, self.num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sampleID_sample_map[i] = cur_sample;
            
    def __shingleThisSampleAndBuildLSHTable(self):
        # Creates a dictionary of sample id - > [list of shingle ids]

        #For now, each sample can only have one shingling.  If the need arises this can be changed.
        
        #sample id to shingle id list
        
        itervar = 0;
                
        # O(M)
        for i in range(0, self.num_of_samples):
            my_sample = self.getSamplebySampleID(i);
            shingles = utils.shingle_sequence(my_sample, self.shingle_size);
            
            shingle_id_for_doc_set = set();
            #O(N)
            for j in range(0, len(shingles)):
                my_shingle = shingles[j];
                
                #Create / Fetch unique shingle ID
                
                #If we've already seen the shingle, retrieve its id for the list
                shingle_id = self.getShingleIDbyShingle(my_shingle);

                if (shingle_id != -1):
                    shingle_id = self.getShingleIDbyShingle(my_shingle);
                    shingle_id_for_doc_set.add(shingle_id);
                    
                #Otherwise, create a new id for it and add it to the known shingle list.
                else:
                    shingle_id = itervar;
                    self.shingleID_shingle_map[shingle_id] = my_shingle;
                    shingle_id_for_doc_set.add(shingle_id);
                
                    itervar = itervar + 1;
                    
                    
                # Now update the shingleID to sampleId dictionary
                k = self.shingleID_to_sampleID_dict.get(shingle_id);
                
                #Do we already have this shingle ID stored from a previous iteration?  If so, add the document to the dictionary - taking care not to duplicate the doc id
                if k:
                    ent = self.shingleID_to_sampleID_dict.get(shingle_id);
                    if (ent.count(i) == 0):
                        ent.append(i);
                        new_entry = [(shingle_id, ent)];   
                        self.shingleID_to_sampleID_dict.update(new_entry);
                
                #otherwise create a new entry.
                else:
                    new_entry = [(shingle_id, [i])];
                    self.shingleID_to_sampleID_dict.update(new_entry);


            shingle_id_for_doc_list = list(shingle_id_for_doc_set);

            update = [(i, shingle_id_for_doc_list)];
            self.sampleID_to_shingleID_dict.update(update);
            self.shingleID_to_sampleID_dict.pop(-1, -1);

    
    def getSampleIdSampleMap(self):
        return self.sampleID_sample_map;
    
    def getLengthDNAString(self):
        return self.length_string;
    
    def getNumOfDNASamples(self):
        return self.num_of_samples;
    
    def getSeed(self):
        return self.seed;
    
    def getAllShingleIds(self):
        return self.shingleID_to_sampleID_dict.keys();
    
    def getAllSamples(self):
                
        sample_list = [];
        for i in range(0, self.num_of_samples):
            sample = self.sampleID_sample_map[i];
            sample_list.append(sample);
        return sample_list;
        
        
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
        return self.sampleID_to_shingleID_dict[sample_id];

    def getShinglesBySampleId(self, sample_id):
        shingle_id_list = self.sampleID_to_shingleID_dict[sample_id];
        shingles_list = [];
        for j in range(0, len(shingle_id_list)):
            my_shingle_id = shingle_id_list[j];
            my_shingle = self.shingleID_shingle_map[my_shingle_id];
            shingles_list.append(my_shingle);
        return shingles_list;
    
    def getShinglesBySample(self, sample_str):
        sample_id = self.getSampleIDbySample(sample_str);
        return self.getShinglesBySampleId(sample_id);

            
    def getShingleSize(self):
        return self.shingle_size;
    
    
    def getSamplesByShingle(self, shingle_str):
        shingle_id = self.getShingleIDbyShingle(shingle_str);

        if (shingle_id == -1):
            return None;
        return self.getSamplesByShingleId(shingle_id);
    
    def getSamplesByShingleId(self, shingle_id):
        sampleIds = self.getSampleIdsByShingleId(shingle_id);
        
        samples = [];
        
        if (sampleIds != None):
            for i in range(0, len(sampleIds)):
                my_sample_id = sampleIds[i];
                sample_str = self.getSamplebySampleID(my_sample_id);
                samples.append(sample_str);
        return samples;
            
    
    def getSampleIdsByShingleId(self, shingle_id):
        sampleIds = self.shingleID_to_sampleID_dict.get(shingle_id);
        return sampleIds;


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
            
        
        
            
        