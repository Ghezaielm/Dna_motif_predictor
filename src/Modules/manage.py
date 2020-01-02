import os
import numpy as np
from datetime import datetime

class Manage():

    def __init__(self):

        # Define dirs paths
        self.curr_dir = os.getcwd()
        self.src_dir = os.path.dirname(self.curr_dir)
        self.app_dir = os.path.dirname(self.src_dir)
        self.data_dir = os.path.join(self.app_dir,os.listdir(self.app_dir)[0])
        self.input_dir = os.path.join(self.data_dir,os.listdir(self.data_dir)[0])
        self.output_dir  = os.path.join(self.data_dir,os.listdir(self.data_dir)[2])
        self.train_input_dir = os.path.join(self.input_dir,os.listdir(self.input_dir)[1])
        self.test_input_dir = os.path.join(self.input_dir,os.listdir(self.input_dir)[0])

        # Initialize models path
        self.model_path = os.path.join(self.data_dir,"models")
        self.model = 0

        # Initialize menu answer
        self.answer = 0

        # is a motif selected ?
        self.state = 0

        # Selected motifs
        # The multi output model is not yet implemented
        # Please make sure to select only one motif
        self.motif_list = []
        self.motif_names = 0
        self.loaded = 0

        # Session attributes
        self.experiment_name = 0

    def menu(self):
        def get_motif_list():
            motif_list = [i[:-4] for i in os.listdir(self.train_input_dir)]
            print("##########")
            print("Motif list")
            print("##########\n\n")
            for id, name in enumerate(motif_list):
                print("{} <- {}".format(id,name))
                self.motif_list.append("{}".format(name))
            return motif_list
        def ask_for_input(motif_list):
            answer = input("Which motif to load ? Unique motif (n), Multiple (n1,n2..]), All (a), Back to menu (q): \n")
            if isinstance(answer,str):
                if answer in ['q','a',"Q",'A']:
                    return answer
                else:
                    try:
                        answer = [int(i) for i in answer.split(",")]
                        self.motif_names = [motif_list[i] for i in answer]
                        self.loaded = 1
                        return answer

                    except Exception as err:
                        print("Invalid entry.")
                        ask_for_input(motif_list)
        def welcome():
            c1 = "1 - Make predictions \n"
            c2 = "2 - Load a model and make predictions\n"
            c3 = "3 - Quit \n"
            if self.loaded==0:
                loaded = ""
                c0 = "0 - Choose your motif(s) {} \n".format(loaded)

            else:
                loaded = "-> motif(s) loaded !"
                c0 = "0 - Choose your motif(s) {} \n".format(loaded)
            print("########################")
            print('DNA Motif Predictor: V1')
            print("########################")
            print("\n","\n",c0,c1,c2,c3)
            ans = input("\n")
            if int(ans)==0:
                self.answer = ask_for_input(get_motif_list())
                welcome()
            if int(ans)==1:
                if self.loaded==0:
                    print("Error, no motif selected\n Choose item 2 to load a previously trained model\n")
                    print(welcome())
                self.state = 1
                self.loaded = 0
                print("Will train a model for:")
                for i in self.motif_names:
                    print("- {}\n\n".format(i))
                now = datetime.now().strftime("%d_%m_%y@%H_%M_%S")

                self.curr_output_dir = os.path.join(self.output_dir, now)
                os.mkdir(self.curr_output_dir)
                self.experiment_name = now
            if int(ans)==2:
                motif_list = np.unique([i.split(".")[0] for i in os.listdir(self.model_path)])
                if len(motif_list)==0:
                    print("Empty")
                    _ = input("Back to menu ? ")
                    welcome()
                else:
                    for id, value in enumerate(motif_list):
                        print(id,"<-",value)

                    model_id = [(i[0],i[1]) for i in enumerate(motif_list)]
                    print("Please, select a model to load")
                    a = input("")
                    if int(a) in [i for i in range(len(model_id))]:
                        self.model_path = os.path.join(self.model_path,model_id[int(a)][1])
                        self.state=2
                        self.motif_names = [model_id[int(a)][1]]
                    else:
                        print("Invalid entry")
                        welcome()

            if int(ans)==3:
                os._exit(1)
            elif int(ans) not in [0,1,2]:
                print("Error, make a new selection.")
                welcome()


        welcome()

