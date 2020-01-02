from data_loader import *

# Data processing
import numpy as np
from datetime import datetime

from sklearn.model_selection import train_test_split
from keras.preprocessing import sequence as sc
from keras.models import model_from_json

class DataProcessing(DataLoader):
    def __init__(self):
        super().__init__()

    @timer
    # Data processing
    def process_data(self):
        '''
         Pad the training data and divide it into training and test set
        :return: pass the processed data as an object attribute
        '''
        print("2) Processing data\n")

        print("Slicing the dataset into training and test sets ...")
        # We aim to train the model with a fraction of the data, or all (default). Throw an error if the fraction is too small or exceed the sample size
        if self.n_inputs > 5:
            # Make sure that the values are sufficient to build a dataset
            try:
                input_data = self.samples[:self.n_inputs]
                input_labels = self.labels[:self.n_inputs]
            except:
                try:
                    input_data = self.samples
                    input_labels = self.labels
                except:
                    print(
                        "Format error: n_input must not exceed sample size.\n Try changing it by using the set_n_input() method.")
        else:
            print(
                "Format error: n_input is too small for training.\n Try changing it by using the set_n_input() method.")

        # Split the data into training and test sets
        X_train, X_test, y_train, y_test = train_test_split(self.samples, self.labels, test_size=self.test_size,
                                                            random_state=42)

        print("Padding the input such as the training sequences have all the same length ...")
        # Store unprocessed training and test sets
        self.X_unp['X_train_unprocessed'] = X_train
        self.X_unp['X_test_unprocessed'] = X_train
        # Truncating and padding the inputs:
        self.X_pro['X_train_processed'] = sc.pad_sequences(X_train, maxlen=self.l_inputs)
        self.X_pro['X_test_processed'] = sc.pad_sequences(X_test, maxlen=self.l_inputs)

        # Store the processed training and test tests
        self.y['y_train'] = np.array([[i] for i in y_train])
        self.y['y_test'] = np.array([[i] for i in y_test])

        del X_train
        del X_test
        del y_train
        del y_test
        del self.samples
        del self.labels

    @timer
    def generate_in_sillico_sequence(self, sequence_length):
        '''
        A function to generate fake input sequences with motifs at specific positions
        :param sequence_length:
        :return: do not return nnothing
        '''
        print(">>> In silico experiment: \n")
        print("Generating a random DNA sequence containing the trained DNA motif.\nLength = {}".format(ls))

        ## We will generate a sequence containing partial and complete motifs
        # First we generate a random sequence of fixed length
        seq = [np.random.choice(self.base, 1)[0] for nt in range(ls)]

        # We will insert the complete motif in 3 places:
        # begining of the seq, middle of the seq, end of the seq
        # We will also insert the partial motif between the begining and the mid of the seq

        ### Create new motifs that contains the consensus
        # First We create a list to store the motifs
        motifs = []
        # store the complete motifs
        for complete in range(3):
            motifs.append([np.random.choice(self.positions[position], 1)[0] for position in self.positions])
        # and the partial ones
        for i in range(2):
            motifs.append([np.random.choice(self.positions[position], 1)[0] for position in self.positions if
                           position < len(self.positions) / 2])
        # convert the motif to numeric sequences
        motifs = [super().convert(seq=sequence) for sequence in motifs]
        # Create a dic to store the motifs sequence positions
        self.exp_dic = {sequence_index: 0 for sequence_index in range(len(motifs))}

        ### Now, we create the fake sequence
        # we start with a random sequence
        seq = [np.random.choice(self.base, 1)[0] for nucleotide in range(sequence_length)]
        # And we construct the TFBDS containing sequence by iterative concatenation of TFBS and random sequences
        for new_seq in range(len(motifs)):
            self.exp_dic[new_seq] = (len(seq), motifs[new_seq])
            seq.extend(motifs[new_seq])
            seq.extend([np.random.choice(self.base, 1)[0] for pos in range(sequence_length)])
        # Finally we convert the input sequence to its reverse complementary (we will also make predictions on that strand)
        self.p_seq = seq
        self.q_seq = super().to_reverse_complementary(seq)
        del seq
        del motifs

    @timer
    def process_inputs(self, seq='fasta'):
        '''
        This function allow to pad the input data so that we can predict motif positions
        :param seq: default is fasta (have to call seq_from_fasta before), else one can provide custom sequence
        :return: pass the input data as attribute
        '''

        print("5) Input processing\n")
        print("Processing input data ...")
        if seq == 'fasta':
            pass
        else:
            self.p_seq = seq
        # First, we store the sequence as an input and convert it to its reverse complementary
        seq = self.p_seq
        self.q_seq = super().to_reverse_complementary(seq)
        # self.l_inputs = len(self.p_seq)
        if self.state==2:
            with open(self.model_path+".txt") as motif_param:
                it = 0
                for line in motif_param:
                    if it == 0:
                        self.l_inputs = [int(i) for i in line.split()][0]
                        it+=1
        # The exp inputs will be composed of a moving windows along the exp sequences

        self.p_inputs = [self.p_seq[i:i + self.l_inputs] for i in range(len(self.p_seq) - self.l_inputs)]
        self.p_inputs = sc.pad_sequences(self.p_inputs, maxlen=self.l_inputs)
        # We do the same on the reverse complementary strand
        self.q_inputs = [self.q_seq[i:i + self.l_inputs] for i in range(len(self.q_seq) - self.l_inputs)]
        self.q_inputs = sc.pad_sequences(self.q_inputs, maxlen=self.l_inputs)

    @timer
    def make_predictions(self, load_model=False):
        '''
         This function allows to make prediction using the trained model or a loaded one
        :param load_model: wheter to load a previously trained model (default is false)
        :return:
        '''
        # We also implemented a function of this function that allows to load the model from the dsik
        print("6) Predictions \n")
        if load_model == False:
            print("Predicting motif positions ...")
            # We make predictions on the direct strand
            self.p_preds = self.model.predict(self.p_inputs)
            # And the same thing in the reverse compelmentary region
            self.q_preds = self.model.predict(self.q_inputs)
            print(self.p_pres, self.q_preds)
        else:

            print("Predicting motif positions ...")
            json_file = open('{}.json'.format(self.model_path), 'r')
            loaded_model_json = json_file.read()
            json_file.close()
            self.model = model_from_json(loaded_model_json)
            # load weights into new model
            self.model.load_weights("{}.h5".format(self.model_path))
            self.p_preds = self.model.predict(self.p_inputs)
            self.q_preds = self.model.predict(self.q_inputs)

    def seq_from_fasta(self, path):
        '''
         This function allows to fetch genomic sequences from a genbank fasta file
        :param path: The filepath
        :return: pass the input fasta sequence as an object attribute
        '''
        if self.state==2:
            now = datetime.now().strftime("%d_%m_%y@%H_%M")
            self.curr_output_dir = os.path.join(self.output_dir, now+"_from_loaded_model")
            print(self.output_dir)
            os.mkdir(self.curr_output_dir)
        self.prediction_dir = os.path.join(self.curr_output_dir,"predictions")
        os.mkdir(self.prediction_dir)
        self.p_seq = []
        with open(path) as fp:
            it = 0
            for line in fp:
                if it >= 1:
                    self.p_seq += super().convert([i for i in line if i != "\n"])
                it += 1

    def free_cache():
        pass