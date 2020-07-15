

$$$$$$$\  $$\   $$\  $$$$$$\        $$\      $$\  $$$$$$\ $$$$$$$$\ $$$$$$\ $$$$$$$$\       
$$  __$$\ $$$\  $$ |$$  __$$\       $$$\    $$$ |$$  __$$\\__$$  __|\_$$  _|$$  _____|      
$$ |  $$ |$$$$\ $$ |$$ /  $$ |      $$$$\  $$$$ |$$ /  $$ |  $$ |     $$ |  $$ |            
$$ |  $$ |$$ $$\$$ |$$$$$$$$ |      $$\$$\$$ $$ |$$ |  $$ |  $$ |     $$ |  $$$$$\          
$$ |  $$ |$$ \$$$$ |$$  __$$ |      $$ \$$$  $$ |$$ |  $$ |  $$ |     $$ |  $$  __|         
$$ |  $$ |$$ |\$$$ |$$ |  $$ |      $$ |\$  /$$ |$$ |  $$ |  $$ |     $$ |  $$ |            
$$$$$$$  |$$ | \$$ |$$ |  $$ |      $$ | \_/ $$ | $$$$$$  |  $$ |   $$$$$$\ $$ |            
\_______/ \__|  \__|\__|  \__|      \__|     \__| \______/   \__|   \______|\__|            
                                                                                            
                                                                                            
                                                                                            
$$$$$$$\  $$$$$$$\  $$$$$$$$\ $$$$$$$\  $$$$$$\  $$$$$$\ $$$$$$$$\  $$$$$$\  $$$$$$$\       
$$  __$$\ $$  __$$\ $$  _____|$$  __$$\ \_$$  _|$$  __$$\\__$$  __|$$  __$$\ $$  __$$\      
$$ |  $$ |$$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ /  \__|  $$ |   $$ /  $$ |$$ |  $$ |     
$$$$$$$  |$$$$$$$  |$$$$$\    $$ |  $$ |  $$ |  $$ |        $$ |   $$ |  $$ |$$$$$$$  |     
$$  ____/ $$  __$$< $$  __|   $$ |  $$ |  $$ |  $$ |        $$ |   $$ |  $$ |$$  __$$<      
$$ |      $$ |  $$ |$$ |      $$ |  $$ |  $$ |  $$ |  $$\   $$ |   $$ |  $$ |$$ |  $$ |     
$$ |      $$ |  $$ |$$$$$$$$\ $$$$$$$  |$$$$$$\ \$$$$$$  |  $$ |    $$$$$$  |$$ |  $$ |     
\__|      \__|  \__|\________|\_______/ \______| \______/   \__|    \______/ \__|  \__|     

Author: GHEZAIEL Morad, Github : Ghezaielm, Mail : ghezaiel.morad@gmail.com 

Here is the first version of DNA Motif Predictor, a python library for DNA motif searching.
 
/!\ Please read the Requirements document to make sure necessary dependencies are well installed

# 1) Introduction and basic use: 

The library is provided with a bunch of training data fetched from the open access JASPAR transcription factor database
and consisting in 746 transcription factor binding motifs from the vertebrata taxa. 
Most of them are given experimental evidences, mainly from Homo Sapiens and Mus Musculus ; please check their website: www.jaspar.genereg.net/
This library allows to train a reccurent neural network (LSTM) to recognize a DNA motif in a given sequence: 

	## a) select a DNA motif from the built-in TF list
	## b) train the model 
	## c) load the DNA sequence for which you want to predict TF binding positions 
	## d) launch prediction 
 
/!\ The multiple TF predictions model is not yet implemented, please make sure to select only one TF. 

You can either use the built-in menu to launch the default pipeline or use the provided class methods throught custom scripts.

# 2) Results 

When training a model for a TF motif, a results folder is created at /data/output and is named as following : day_month_year@hour_minutes_second
It will contain the model error metric during training and a subfolder named predictions.
Each models are automatically saved in the /data/models folder. 
When loading a previously trained module, your predictions are stored as results from a new experiment. 
You can load your sequence to predict in fasta format and store it in /data/input/test.

# 3) Important 

The library is provided with a default input sequence, Nanog human (Genbank). 
To make new predictions (ie with new sequences) throught the default menu, please make sure to: 

- store the fasta sequence file in the /data/input/test folder 
- replace the path accordingly in the main script 
- Add the annotations (start,end,intron) manually in the dictionnary (see main script) 



 







                                                                                            
                                                                                            
