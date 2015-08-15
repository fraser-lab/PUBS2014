# Needs the following pickles in the pkl folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Needs pkls to be analyzed in fitness folder.

Module 3: Outlier detection and removal
toss_outliers.py Barcode_scores.pkl codon 4 -o clean_BCs.pkl

clean_BCs.pkl
	{“dirty_barcodes”: [(aa_position_in_Ub, Mutant_Codon, barcode) …]
	“clean_barcodes”: [(aa_position_in_Ub, Mutant_Codon, barcode) …]}
	“dirty_barcodes”: [(74, AUU, “AGCTCTA”), (21, CCC, “ACTTCTA”) …]
	“clean_barcodes”: [(5, GCU, “AGCCCTA”), (21, UUU, “GCATTTC”) …]}

This script compares the scores of barcodes mapping to the same allele. The median absolute deviation (MAD) is calculated for the barcodes that map to the same allele. Outlier barcodes are determined as scores that are greater than or equal to 1.5 times the interquartile range of the distribution and removed. The codon flag tells the script to compare all scores mapping to the same codon. The 4 flag sets the minimum number of barcodes before the MAD will be performed. The output is a pkl of barcodes to be kept in the dataset.


remove_bad_BCs.py Barcode_scores.pkl clean_BCs.pkl Allele_Dictionary.pkl -o Barcode_scores_outliers_removed.pkl 
This script checks the Barcode_scores dictionary against the clean_BCs returned by toss_outliers.py and removes those BCs that are not in the clean_BCs.pkl. The script returns a dictionary keyed on sample barcode with fitness scores as values but with outlier BCs removed.

Barcode_scores_outliers_removed.pkl 
	{ “barcode” : float(Fitness_score)}
	“AGCCCTA” : -0.1
	“GCATTTC” : -0.67}
