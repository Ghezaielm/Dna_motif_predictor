from data_processing import *
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers.embeddings import Embedding
import matplotlib.pyplot as plt

class Model(DataProcessing):
    def __init__(self):
        # The model
        super().__init__()
        self.model = "unknown"
        # The model history
        self.history = "unknown"
        # Loss
        self.loss = 'binary_crossentropy'
        # Activation
        self.activation = 'sigmoid'
        # Optimizer
        self.optimizer = 'adam'
        # Metrics
        self.metrics = ['accuracy']
        # LSTM units
        self.n_LSTM_units = 50
        # Encoding space
        self.encoding_space_size = 64
        # Training steps
        self.epochs = 100
        # Batchsize
        self.batch_size = 70

    def set_test_size(self, value):
        # We set the testing size
        self.test_size = value

    def set_loss(self, loss):
        # We set the loss
        print("Warning: alternative losses are not yet implemented, make sure the loss is setted to 'adam'.")
        self.loss = loss

    def set_optimizer(self, opt):
        # We set the optimizers
        self.optimizer = opt

    def set_metrics(self, metrics):
        # We set the metrics
        self.metrics = metrics

    def set_LSTM_units(self, n):
        # We set the lstm units (number)
        self.n_LSTM_units = n

    def set_encoding_size(self, n):
        # And finally the encoding size ( ie the size of the embedding space)
        self.encoding_space_size = n

        # Instanciate the model with the params

    def set_model(self):
        '''
        This function allows to set the model
        :return: pass the model as an object attribute
        '''
        # Compiling the model with such loss
        evl = self.encoding_space_size
        # And we create the model
        self.model = Sequential()
        # We set the embedding layer : convert the input into a new space
        self.model.add(Embedding(self.n_inputs, evl, input_length=self.l_inputs))
        # We add LSTMs units
        self.model.add(LSTM(self.n_LSTM_units))
        # And finally, a dense layer with a sigmoid activation to output probabilities
        self.model.add(Dense(1, activation=self.activation))
        # Now we compile the model with the params
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)
        # We print the model parameters
        print(self.model.summary())
        # Store the model as an attribute of the experiment object

    @timer
    # We train the model with the data
    def train_model(self):
        '''
         This function allows to train the model
        :return: pass the trained model as an object attribute
        '''
        print("3) Training the model\n")
        print(input("Ready ?:\n"))
        from keras.callbacks import EarlyStopping
        print(
            ">>> Training: \nTraining size: {}\nTesting size: {}\nTrain for: {} epochs.\nSamples by update: {}\nMetric: {}".format(
                int(len(self.X_pro['X_train_processed']) * (1 - self.test_size)),
                int(len(self.X_pro['X_train_processed']) * self.test_size), self.epochs, self.batch_size, self.metrics))
        # And fit the model with the data
        es = EarlyStopping(monitor='accuracy', mode='max',restore_best_weights=True,patience=5, verbose=1)
        self.history= self.model.fit(self.X_pro['X_train_processed'], self.y['y_train'], epochs=self.epochs,
                                 batch_size=self.batch_size, callbacks=[es],verbose=0)

        # And we plot the value
        loss = self.history.history["loss"]
        plt.style.use('dark_background')
        f = plt.figure(figsize=(10, 7))
        ax = f.add_subplot(111)
        ax.set_title("Model error during training", fontsize=20)
        ax.set_xlabel("epochs", fontsize=20)
        ax.set_ylabel("Error", fontsize=20)
        ax.plot(loss, "-o", linewidth=3)
        if self.state == 0:
            now = datetime.now().strftime("%d_%m_%y@%H_%M")
            self.curr_output_dir = os.path.join(self.output_dir, now)
            os.mkdir(self.curr_output_dir)

        plt.savefig(os.path.join(self.curr_output_dir,"model_error.png"))
        plt.show()
        # Finally, we store the model as an attribute
        model_json = self.model.to_json()
        self.model_path = os.path.join(self.model_path,self.motif_names[0])
        with open("{}.json".format(self.model_path), "w") as json_file:
            json_file.write(model_json)
        # Serialize weights to HDF5
        self.model.save_weights("{}.h5".format(self.model_path))
        print("Model saved to disk")

    # We make predictions on the test set
    def evaluate_model(self):
        '''
         This function allows to evaluate the model by launching prediction on the test set
        :return: output the metrics in stdout
        '''
        y_preds = self.model.predict(self.X_pro['X_test_processed'])
        y_true = [i[0] for i in self.y['y_test']]
        classes = [1 if i < 0.8 else 0 for i in y_preds]
        map = [i for i in zip(classes, y_true)]
        VP = sum([1 for i in y_true if i == 1])
        FP = sum([1 for i in map if i == (1, 0)])
        FN = sum([1 for i in map if i == (0, 1)])
        def precision():
            return VP/(VP+FP)

        def recall():
            return VP / (VP + FN)

        def f1():
            return 2*(precision()*recall())/(precision()+recall())

        # We evaluate how well the model perform on unseen data
        scores = self.model.evaluate(self.X_pro['X_test_processed'], self.y['y_test'], verbose=0)
        print("4) Model evaluation\n")
        print("Performing cross validation ...")
        # And print the metrics
        print("Evaluation accuracy: {}".format((scores[1] * 100)))
        print("Precision: {}".format(round(precision()*100,2)))
        print("Recall: {}".format(round(recall()*100,2)))
        print("F1-score: {}".format(round(f1()*100,2)))
        print("\n")