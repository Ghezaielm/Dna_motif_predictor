from manage import *

# Data processing
from multiprocessing import Pool
import numpy as np
import time
import os

# PEP 257: Docstrings

def timer(func):
    '''

    :param func: a function to call for which we want to measure the ecution time
    :return: a wrapper


    '''
    def wrapper(*args):
        start = time.time()
        func(*args)
        end = time.time() - start
        print("Time elapsed: {} sec(s)\n".format(round(end, 2)))

    return wrapper


class DataLoader(Manage):
    # This class allows to train a RNN model to predict the presence of a motif within a DNA sequence
    def __init__(self):

        ######################################################
        #######          General attributes             ######
        ######################################################
        # Inherite session attributes
        super().__init__()
        #self.experiment_name = "Unknown"
        self.creation_date = "Unknown"
        self.mode = 0

        ######################################################
        #####          Training data attributes          #####
        ######################################################

        ## Nucleotides frequencies for TF binding sites (Default is SMAD Human)
        self.A = [166, 196, 14, 33, 0, 0, 106, 211, 0, 765, 30, 117, 110]
        self.T = [432, 75, 8, 45, 899, 0, 231, 125, 896, 19, 736, 583, 214]
        self.G = [168, 34, 877, 103, 0, 0, 539, 280, 3, 18, 133, 47, 35]
        self.C = [133, 594, 0, 718, 0, 899, 23, 283, 0, 97, 0, 152, 540]
        ## Normalized nt frequency matrix
        self.normalized_matrix = {i:[] for i in range(4)}
        ## Number of training inputs
        self.n_inputs = "Unknown"
        ## Motif length
        self.l_inputs = "Unknown"
        ## The encoding base
        self.base = [1, 3, 5, 7]
        # Training sequences
        self.samples = "Unknown"
        self.null_samples = []
        # The training labels
        self.labels = "Unknown"
        self.null_labels = []
        # A map of nucleotide frequencies at each motif position
        self.positions = "Unknown"
        # Fraction of the training data to use as test set
        self.test_size = 0.33
        # Copy of the data in dics
        self.X_unp = {'X_train_unprocessed': 0, 'X_test_unprocessed': 0}
        self.X_pro = {'X_train_processed': 0, 'X_test_processed': 0}
        self.y = {'y_train': 0, "y_test": 0}

        ##########################################################
        #######            Input data attributes          ########
        ##########################################################

        # a dic containing nt freq for the experimental set
        self.exp_dic = "Unknown"
        # The experimental sequence (direct strand)
        self.p_seq = "Unknown"
        # The experimental sequence (reverse strand)
        self.q_seq = "Unknown"
        # inputs for the direct strand
        self.p_inputs = "Unknown"
        # inputs for the reverse strand
        self.q_inputs = "Unknown"
        # Predictions for the direct strand
        self.p_preds = "Unknown"
        # Predictions for the reverse strand
        self.q_preds = "Unknown"


    ###########################################################
    #####                     Setters                     #####
    ###########################################################

    def set_motif_frequencies(self, matrix):
        '''
        Load the experiment object with a custom motif
        :param matrix:a matrix of nucleotide frequencies in format base * position
        :return:
        '''
        # This function allows to set the nucleotides frequencies in order to build training sequences
        self.A = matrix[0]
        self.T = matrix[1]
        self.G = matrix[2]
        self.C = matrix[3]

    def set_n_inputs(self, n):
        # We set the n_inputs parameters as attributes
        self.n_inputs = n

    def set_base(self, base):
        # We set the encoding base
        self.base = base

    # Training the model with chipseq data (to implement)
    # Set sequences
    def set_samples(self, samples):
        self.samples = samples

    # Set labels (bound/not bound)
    def set_labels(self, labels):
        if np.unique(labels) == [0]:
            self.labels = labels
        else:
            print("Invalid labels, make sure they are binary encoded")
    # Set the test size
    def set_test_size(self, value):
        if value < 1 and value > 0:
            self.test_size = value
        else:
            print("Invalid test size, try to set it to 0.33")

    def set_motif_names(self,name):
        if type(name)==str:
            self.motif_names = [name]
        else:
            print('Invalid motif name, launch the menu to get the motif list or browse the /input/train folder')

    ##############################################################
    #####                       Getters                      #####
    ##############################################################

    # Get the motif(s) name(s)
    def get_motif_names(self):
        ''' Return motif names'''
        [print(i) for i in self.motif_list]

    # Get the nt frequencies matrix
    def get_matrix(self,mode):
        ''' Get the motif nucleotide frequency matrix

        :param mode: 'unique' (str) # multiple motif mode is not yet implemented
        :return: matrix of nucleotide frequency passed as an attribute
        '''
        if mode=="unique":
            with open(os.path.join(self.train_input_dir,"{}.txt".format(self.motif_names[0]))) as f:
                dic = {"A":[],"T":[],"G":[],"C":[]}
                base_iter = 0

                for line in f:
                    if [i for i in line][0]=='1':
                        pass
                    else:
                        line = line.split(",")
                        if base_iter ==0:
                            curr = 'A'
                        else:
                            curr = line[0]

                        dic[curr]= [int(i) for i in line[1:]]
                        base_iter+=1
                        if base_iter ==4:
                            break

                self.A = dic["A"]
                self.T = dic["T"]
                self.G = dic["G"]
                self.C = dic["C"]

                for i in range(len(self.A)):
                    if self.A[i]==0 and self.T[i]==0 and self.G[i]==0 and self.C[i]==0:
                        self.A[i] = 10
                        self.T[i] = 10
                        self.G[i] = 10
                        self.C[i] = 10


    def get_strict_motif(self):
        '''
        :return: return the consensus motif from the experimental object
        '''
        motif = []
        [motif.append("N") if len(np.unique([self.A[i], self.T[i], self.G[i], self.C[i]])) == 1
         else motif.append({0: "A", 1: "T", 2: "G", 3: "C"}[[self.A[i], self.T[i], self.G[i], self.C[i]].index(
            max([self.A[i], self.T[i], self.G[i], self.C[i]]))])
         for i in range(len(self.A))]
        return motif

    def get_params(self):
        '''

        :return: return the experiment details
        '''
        if self.experiment_name == "Unknown":
            print("Cache is empty: No experiment found.")
        else:
            print("Experiment name {}".format(self.experiment_name))
            print("##### Default params:")
            for attributes, value in self.__dict__items():
                if isinstance(value, (int, str, float)):
                    print("> ", attributes, " = ", value)
                else:
                    print("> ", attributes, " = ", "Empty")


    ############################################################
    #####                   Methods                       ######
    ############################################################

    @timer
    def normalize_matrix(self):
        '''

        :return:
        '''
        for i in range(len(self.A)):
            v = [self.A[i],self.T[i],self.G[i],self.C[i]]
            v = [i/max(v) for i in v]
            [self.normalized_matrix[j].append(v[j]) for j in range(len(v))]



    def get_alignment_score(self,seq):
        if type([seq[0]])==str:
            seq = DataLoader.convert(self,seq)
            ind  = {1:0,3:1,5:2,7:3}
            seq = [ind[i] for i in seq ]
        else:
            ind = {1: 0, 3: 1, 5: 2, 7: 3}
            seq = [ind[i] for i in seq]
        scores = [self.normalized_matrix[seq[i]][i] for i in range(len(seq)) ]
        return sum(scores)/len(seq)

    def count_to_frequency(self):
        '''
        Convert the count matrix of the experimental object frequencies
        :return: pass the nucleotide frequency matrix as experimental attribute
        '''
        self.l_inputs = len(self.A)
        with open("{}.txt".format(os.path.join(self.model_path,self.motif_names[0])), "w") as motif_param:
            motif_param.write(str(self.l_inputs))
        # We create a dic to store the positions
        positions = {i: [] for i in range(len(self.A))}
        # For each nucleotides positions
        for position in range(len(self.A)):
            # We add the data to a list
            nucleotides = [self.A, self.T, self.G, self.C]
            # We store the nts frequencies inside a list
            curr = [nucleotide[position] if nucleotide[position] > 0 else 0 for nucleotide in nucleotides]
            # We convert theses frequencies to probabilities
            curr = [round(freq / sum(curr), 1) if freq != 0 else 0 for freq in curr]

            # We create here a list containing nucleotides at various frequencies for random picking

            for nucleotide in range(len(curr)):
                # for each nucletides frequency
                val = 100 * curr[nucleotide]
                # We add k times the corresponding nucleotide
                for k in range(int(val)):
                    if nucleotide == 0 and val > 0:
                        positions[position].append("A")
                    if nucleotide == 1 and val > 0:
                        positions[position].append("T")
                    if nucleotide == 2 and val > 0:
                        positions[position].append("G")
                    if nucleotide == 3 and val > 0:
                        positions[position].append("C")
                # Resulting in position specific list that contains A, T, G and C 'k times' such as its len is equal to 100 (for random picking)

        # We store the nt frequencies for each position as an object attribute
        self.positions = positions

    # NT to base
    def convert(self, seq):
        return [{"A": self.base[0], "T": self.base[1], "G": self.base[2], "C": self.base[3]}[i] for i in seq]

    # Base to NT
    def reverse(self, seq):
        return [{self.base[0]: "A", self.base[1]: "T", self.base[2]: "G", self.base[3]: "C"}[i] for i in seq]

    # base to reverse complementary
    def to_reverse_complementary(self, seq):
        dic = {1: 3, 3: 1, 5: 7, 7: 5}
        return [dic[i] for i in seq]

    @timer
    # Generate random DNA sequences containing the binding site
    def generate_sequences(self, n, thresh_label):
        '''
         This function generates training data by sampling from the frequency matrix
        :param
        :param n: The number of exemples to sample
        :param thresh_label: correlation between positive and negative samples
        :return: pass the training data as attributes
        '''

        print("1) Loading data\n")
        print("Generating training sequences ...")
        # Store the lenght of the generate sequences as an attribute
        self.n_inputs = n
        thresh_label = [-thresh_label,thresh_label]
        # We generate random sequences containing the motif according to nucleotides frequencies
        self.samples = [[np.random.choice(self.positions[position], 1)[0]
                        for position in self.positions]
                        for new_input in range(self.n_inputs)]


        #### Generate sequences that does not contain the motif (negative samples)
        # Convert the dataset in a format handl able by the model
        # Let's assign a positive label to our sequences:
        self.labels = [1 for sample in self.samples]
        print("Creating random sequences that do not contain the DNA motif ...")
        # Now, lets create a set of sequences with null labels


        try:
            self.samples = [DataLoader.convert(self, sample) for sample in self.samples]
        except Exception:
            pass
        # error point

        def generate_null_samples(threaded):

            for sample in threaded:
                # For each sample
                state = 0
                # While we still have not generated a sequence that is not correlated to the computed sample
                while state == 0:
                    # We generate a random sequence
                    sequence = [np.random.choice(sample, 1)[0] for nt in sample]
                    # And make sure that its not correlated to the sample
                    corr = np.corrcoef(sequence, sample)
                    # If the we find an uncorrelated sequence (ie: does not contain the motif)
                    if corr[0][1] > thresh_label[0] and corr[0][1] < thresh_label[1]:
                        # We keep that data to construct a balanced dataset (negative exemple)
                        self.null_samples.append(sequence)
                        self.null_labels.append(0)
                        # We update the state
                        state += 1

        #############################################################################################################
        # One can call the above function in a parallelized way, not yet implemented but here are the bribes
        try:
            with Pool(5) as pool:
                pool.map(generate_null_samples,self.samples)

        except Exception:
            generate_null_samples(self.samples)
        """    
        self.thread_size = 4
        to_thread = [self.samples[int(len(self.samples)* i):int(len(self.samples)* (i+1/self.thread_size))] for i in np.arange(0,1,1/self.thread_size)]
        threads = [Thread(target=generate_null_samples, args=(i,)) for i in to_thread]
        [i.start() for i in threads]
        [i.join() for i in threads]"""
        ################################################################################################################

        # We add the resulting data to the core dataset
        self.samples += self.null_samples
        self.labels += self.null_labels
        del self.null_samples
        del self.null_labels
