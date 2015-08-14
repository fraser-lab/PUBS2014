Module 2: Initial scoring - Barcodes, initial counts cutoff = 3
Script1:
pickle_condensing.py TS_1.pkl TS_2.pkl TS_3.pkl perturbation replicate -o Barcode_Counts.pkl

This script simply takes the counts from the dictionaries created by Module 0 Script 3 and combines them into a single dictionary that contains the counts for a given barcode for all 3 samples that describe an experiment. The perturbation and replicate inputs are used in naming the output dictionary.

Script2:
Score_BCs.py Barcode_Counts.pkl time1 time2 -o Barcode_scores.pkl 

Barcode_Counts.pkl
	{ “barcode” : [count_sample1, count_sample2…]}
	“AGCTCTA” : [15, 3, 1]
“AGCCCTA” : [222, 23, 21]}


This script takes a .pkl of a dictionary (Barcode_Counts.pkl) keyed on the sample barcodes with values of a list of counts at each time point. The scoring function uses these counts and scores them based on the time values (in WT generations) - the relative fitness is compared to wild type barcodes, which are distinguished in Allele_Dictionary.pkl. Output is a dictionary keyed on sample barcode with values of the fitness score. Fitness scores are determined by calculating the slope of the regression line of the three counts for each barcode. The score reported is log2(Mutant Slope/ WT Slope). Any barcode that is observed three or less times in the initial sample is not used in the fitness score calculation

Barcode_scores.pkl
	{ “barcode” : float(Fitness_score)}
	“AGCTCTA” : -0.56
	“AGCCCTA” : -0.1}
