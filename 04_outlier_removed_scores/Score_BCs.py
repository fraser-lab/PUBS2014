'''
This script takes a .pkl of a dictionary (Barcode_Counts.pkl) keyed on the sample barcodes with values 
of a list of counts at each time point. The scoring function uses these counts and scores them based on
the time values (in WT generations) - the relative fitness is compared to wild type barcodes, which are 
distinguished in Allele_Dictionary.pkl. Output is a dictionary keyed on sample barcode with values of the 
fitness score. Fitness scores are determined by calculating the slope of the regression line of the three 
counts for each barcode. The score reported is log2(Mutant Slope/ WT Slope). Any barcode that is observed 
three or less times in the initial sample is not used in the fitness score calculation

usage:
Score_BCs.py Barcode_Counts.pkl time1 time2

output = Barcode_scores.pkl 
'''

import cPickle as pic
from scipy import stats
import numpy as np
import collections
import sys
import math
import matplotlib.pyplot as plt
import os

data = pic.load(open(sys.argv[1],'r'))
t1 = float(sys.argv[2])
t2 = float(sys.argv[3])

TRANSLATE = {'ACC': 'T', 'GUC': 'V', 'ACA': 'T', 'AAA': 'K', 'GUU': 'V', 'AAC': 'N', 'CCU': 'P', 'UAU': 'Y', 'AGC': 'S', 'CUU': 'L', 'CAU': 'H', 'AAU': 'N', 'ACU': 'T', 'GUG': 'V', 'CAC': 'H', 'ACG': 'T', 'AGU': 'S', 'CCA': 'P', 'CCG': 'P', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'AUG': 'M', 'UGC': 'C', 'CAG': 'Q', 'UGA': 'STOP', 'UGG': 'W', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': 'STOP', 'GGA': 'G', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'AUC': 'I', 'GGC': 'G', 'GCG': 'A', 'CGC': 'R', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AUU': 'I', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'CGA': 'R', 'UAG': 'STOP', 'GCU': 'A', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}

cutoff = 3

timeList = [0, t1, t2]

T0WTCounts = 0.
T1WTCounts = 0.
T2WTCounts = 0.
for key in data.keys():
	if key[1] == 'WT':
		T0WTCounts += data[key][0]
		T1WTCounts += data[key][1]
		T2WTCounts += data[key][2]
wt_counts = [T0WTCounts,T1WTCounts, T2WTCounts]



ratio_dict = {}
all_zeros = 0
for key in data.keys():
	mut_counts = data[key]
	if mut_counts == [0 for i in range(len(mut_counts))]:
		all_zeros +=1
		pass
	else:
		if mut_counts[0] >= cutoff:
			mut_counts[0]
			ratio_dict[key] = [] 
			for i in range(len(mut_counts)):
				if mut_counts[i] == 0:
					mut_counts[i] = 0.5
				else:
					mut_counts[i] = float(mut_counts[i])
				ratio_dict[key].append(math.log(float(mut_counts[i]/float(wt_counts[i])),2)) 



score_dict = {}

for key in ratio_dict:
	slopeOverall, interceptOverall, r_valOverall, p_valOverall, std_errOverall = stats.linregress(timeList, ratio_dict[key])
	score_dict[key] = slopeOverall

out_name = 'Barcode_scores' + str(sys.argv[1])
pic.dump(score_dict, open(out_name,'w'))
