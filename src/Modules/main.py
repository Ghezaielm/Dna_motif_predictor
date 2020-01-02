from show import *
import time


if __name__ == '__main__':

        # Experiment instanciation
        st = time.time()
        exp = Experiment()
        exp.menu()

        # Whether to load a previous model
        if exp.state==2:
                exp.seq_from_fasta(exp.test_input_dir+"\\nanog.txt")
                exp.process_inputs('fasta')
                exp.make_predictions('True')
                dic = {"Nanog": [214, 5694, 'p', [[i for i in range(364, 3549)],
                                                  [i for i in range(3811, 5056)],
                                                  [i for i in range(5142, 5278)]]]}
                exp.show_positions('fasta', 10, dic)

        # Or launch the basic pipeline
        else:

                # exp.get_params()
                exp.get_matrix(mode="unique")
                exp.count_to_frequency()
                exp.generate_sequences(1000,0.1)
                exp.process_data()
                exp.set_model()
                exp.train_model()
                exp.evaluate_model()
                exp.seq_from_fasta(exp.test_input_dir+"\\nanog.txt")
                exp.process_inputs()
                exp.make_predictions('False')

                dic = {"Nanog": [214, 5694, 'p',[[i for i in range(364,3549)],[i for i in range(3811,5056)],[i for i in range	(5142,5278)]]]}
                exp.show_positions('fasta', 10, dic)

                print("Total Run time: {} sec(s)".format(int(time.time() - st)))

