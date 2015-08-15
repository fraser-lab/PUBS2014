'''
script to take the output of the outlier removeal (a list of clean BCs) and remove them from the dictionaries
output byt the pickle condensing script. Outputs a dictionary with bad BCs removed to send to fitness calculation
'''
import cPickle as pic
from scipy import stats
import numpy as np
import collections
import sys
import math
import matplotlib.pyplot as plt

input_counts = pic.load(open(sys.argv[1], 'rb')) # input dictionary is the origianl pickle condensing output
clean_dic = pic.load(open(sys.argv[2], 'rb')) #output of outlier dectection
allele_dict = pic.load(open('../../pkl/allele_dic_with_WT.pkl', 'r'))
output_counts = {}

for BC in clean_dic['clean_barcodes']:
	output_counts[BC[2]] = input_counts[BC[2]]
	#output_counts[BC] = input_counts[BC]

for BC in input_counts.keys():
	if allele_dict[BC][1] == 'WT':
		output_counts[BC] = input_counts[BC]


pic.dump(output_counts, open((str(sys.argv[1][:-4]) + '_outliers_removed.pkl'), 'wb'))