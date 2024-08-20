import utils;
import TwoWayDict;
import editdistance;
import math;
import random;

class DnaSequenceDataSampleAndLSHTable:


    def __init__(self, length_string, num_of_samples, shingle_size, sampling_choice, num_dna_seeds = -1, num_of_changes_per_seq = -1,  seed=-1):
        #Sampling choice 1 requires no additional params
        #Sampling choice 2 requires maximum_edit_distance, num_dna_seeds, num_changes_per_seq
        #Sampling choice 3 requires maximum_edit_distance
        self.length_string = length_string;
        self.num_of_samples = num_of_samples;
        self.seed = seed;
        self.shingle_size = shingle_size;
        self.num_random_dna_seeds = num_dna_seeds;
        self.sampling_choice = sampling_choice;
        self.num_of_changes_per_seq = num_of_changes_per_seq;

        self.sampleID_sample_map = TwoWayDict.TwoWayDict();
        self.shingleID_shingle_map = TwoWayDict.TwoWayDict();
        self.sampleID_to_shingleID_dict = {};
        #this is the LSH table.  Needs an initial entry to avoid key error.
        self.shingleID_to_sampleID_dict = {-1 : -1 };
        self.all_sample_ids = [];


        #Basic random sample
        if (sampling_choice == 1):
            self.__generateRandomSamples();
        #Guided random sample
        elif(sampling_choice == 2):
            self.__generateGuidedRandomSamples();
        #Centered random sample
        elif(sampling_choice == 3):
            self.__generateCenteredRandomSample();
        else:
            raise Exception("Bad sampling choice passed.  Sampling choice int:  ", sampling_choice, " valid values are 1,2,3")

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
            sample_list = self.getSamplesByShingle(my_shingle);

            if (sample_list != None):
                out_set.update(sample_list);

        out_samples = list(out_set);
        return out_samples;

    def getLocalSampleIds(self, query_dna_sequence):
        #Given a query DNA sequence, return all 'local' samples based on the LSH scheme.
        shingles = utils.shingle_sequence(query_dna_sequence, self.shingle_size);
        out_set = set();

        for i in range(0, len(shingles)):
            my_shingle = shingles[i];
            sample_id_list = self.getSampleIdsByShingle(my_shingle);

            if (sample_id_list != None):
                out_set.update(sample_id_list);

        out_samples = list(out_set);
        return out_samples;


    def __generateRandomSamples(self):
        for i in range(0, self.num_of_samples):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sampleID_sample_map[i] = cur_sample;

    def __generateGuidedRandomSamples(self):
        iter_var = 0;

        #Put the seed DNA samples at the front of the dictionary
        #That's why we do this for loop twice
        for i in range(0, self.num_random_dna_seeds):
            cur_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
            self.sampleID_sample_map[iter_var] = cur_sample;
            iter_var = iter_var + 1;

        self.num_replicas_per_seed = math.floor(self.num_of_samples/self.num_random_dna_seeds);
        remaining_seeds = self.num_of_samples - self.num_random_dna_seeds;

        for i in range(0, self.num_random_dna_seeds):
            cur_sample = self.getSamplebySampleID(i);
            remaining_seeds = self.num_of_samples - self.num_random_dna_seeds;
            amt_seeds_per_changes = math.floor(remaining_seeds / self.num_of_changes_per_seq);
            for j in range(0, amt_seeds_per_changes):
                #Edit distance computation is very expensive, so we won't do it here.
                my_seq = utils.change_2N_edit_distance_in_DNA_seq(cur_sample, self.num_of_changes_per_seq);
                self.sampleID_sample_map[iter_var] = my_seq;
                iter_var = iter_var + 1;

        #If there's any remaining due to the floor operation, just start populating
        while (iter_var < self.num_of_samples):
            samp_inx = random.randint(0, self.num_random_dna_seeds);
            cur_sample = self.getSamplebySampleID(samp_inx);
            my_seq = utils.change_2N_edit_distance_in_DNA_seq(cur_sample, self.num_of_changes_per_seq);
            self.sampleID_sample_map[iter_var] = my_seq;
            iter_var = iter_var + 1;

    def __generateCenteredRandomSample(self):

        #        self.max_num_of_changes_per_seq = num_of_changes_per_seq;

        #Not implemented yet
        centered_sample = utils.generate_random_dna_sequence(self.length_string, self.seed);
        self.sampleID_sample_map[0] = centered_sample;

        iter_var = 1;
        vals_each_distance = math.floor(self.num_of_samples / self.num_of_changes_per_seq);

        for j in range(1, self.num_of_changes_per_seq):
            for r in range(0, vals_each_distance):
                my_seq = utils.change_2N_edit_distance_in_DNA_seq(centered_sample, j);
                self.sampleID_sample_map[iter_var] = my_seq;
                iter_var = iter_var + 1;

        #Fill the rest with random distances
        while (iter_var < self.num_of_samples):
            my_dist = random.randint(1, self.max_num_of_changes_per_seq);
            my_seq = utils.change_2N_edit_distance_in_DNA_seq(centered_sample, my_dist);
            self.sampleID_sample_map[iter_var] = my_seq;
            iter_var = iter_var + 1;


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

    def getOneWaySampleIdSampleMap(self):

        output = {};

        update_list = [];
        for i in range(0, self.num_of_samples):
            update_list.append([i, self.sampleID_sample_map[i]]);

        output.update(update_list);
        return output;

    def getTwoWaySampleIdSampleMap(self):
      return self.sampleID_to_shingleID_dict;

    # This is the number of initial samples used for guided random sampling
    def getNumDnaSeeds(self):
        return self.num_random_dna_seeds;

    def getMaxNumChangesPerSequence(self):
        return self.max_num_of_changes_per_seq;


    def getLengthDNAString(self):
        return self.length_string;

    def getNumOfDNASamples(self):
        return self.num_of_samples;

    def getSeed(self):
        return self.seed;

    def getAllShingleIds(self):
        return list(self.sampleID_to_shingleID_dict.keys());

    def getAllSampleIds(self):
        return list(self.sampleID_to_shingleID_dict.keys());

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

    def getSampleIdsByShingle(self, shingle_str):
        shingle_id = self.getShingleIDbyShingle(shingle_str);

        if (shingle_id == -1):
            return None;
        return self.getSampleIdsByShingleId(shingle_id);




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
