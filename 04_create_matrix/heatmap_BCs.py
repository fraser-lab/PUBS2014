import cPickle as pic
from collections import defaultdict
import pylab
import operator
import sys
import numpy as np
import string
import matplotlib.pyplot as plt

dictionary = pic.load( open( sys.argv[1], "rb" ) )
allele_dict = pic.load(open(sys.argv[2],'rb'))
outfile = sys.argv[3]

fit_matrix = np.zeros((21,76))
dev_matrix = np.zeros((21,76))
table = string.maketrans('T', 'U')
TRANSLATE = {'ACC': 'T', 'GUC': 'V', 'ACA': 'T', 'AAA': 'K', 'GUU': 'V', 'AAC': 'N', 'CCU': 'P', 'UAU': 'Y', 'AGC': 'S', 'CUU': 'L', 'CAU': 'H', 'AAU': 'N', 'ACU': 'T', 'GUG': 'V', 'CAC': 'H', 'ACG': 'T', 'AGU': 'S', 'CCA': 'P', 'CCG': 'P', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'AUG': 'M', 'UGC': 'C', 'CAG': 'Q', 'UGA': 'STOP', 'UGG': 'W', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': 'STOP', 'GGA': 'G', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'AUC': 'I', 'GGC': 'G', 'GCG': 'A', 'CGC': 'R', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AUU': 'I', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'CGA': 'R', 'UAG': 'STOP', 'GCU': 'A', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, 'STOP': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}

def average_list(list):
    return float(sum(list)/len(list))
    
TRANSLATE_num = {"UUU":"1", "UUC":"2", "UUA":"3", "UUG":"4",
    "UCU":"5", "UCC":"6", "UCA":"7", "UCG":"8",
    "UAU":"9", "UAC":"10", "UAA":"11", "UAG":"12",
    "UGU":"13", "UGC":"14", "UGA":"15", "UGG":"16",
    "CUU":"17", "CUC":"18", "CUA":"19", "CUG":"20",
    "CCU":"21", "CCC":"22", "CCA":"23", "CCG":"24",
    "CAU":"25", "CAC":"26", "CAA":"27", "CAG":"28",
    "CGU":"29", "CGC":"30", "CGA":"31", "CGG":"32",
    "AUU":"33", "AUC":"34", "AUA":"35", "AUG":"36",
    "ACU":"37", "ACC":"38", "ACA":"39", "ACG":"40",
    "AAU":"41", "AAC":"42", "AAA":"43", "AAG":"44",
    "AGU":"45", "AGC":"46", "AGA":"47", "AGG":"48",
    "GUU":"49", "GUC":"50", "GUA":"51", "GUG":"52",
    "GCU":"53", "GCC":"54", "GCA":"55", "GCG":"56",
    "GAU":"57", "GAC":"58", "GAA":"59", "GAG":"60",
    "GGU":"61", "GGC":"62", "GGA":"63", "GGG":"64",}


aa_dict = {}
for key in dictionary:
    codon = allele_dict[key]
    if codon[1] != 'WT':
        aa_num = codon[0]
        codon = codon[1]
        codon = string.translate(codon, table)
        aa = str(TRANSLATE[codon])
        aa_val = int(transform[aa])
        fit_score = dictionary[key][0]
        aa_key = (aa_num, aa)
    if aa_key in aa_dict.keys():
        aa_dict[aa_key].append(fit_score)
    else:
        aa_dict[aa_key] = [fit_score]

for aa in aa_dict.keys():
    aa_dict[aa] =  [average_list(aa_dict[aa]), np.std(aa_dict[aa])]


for aa in aa_dict.keys():
    fitness = aa_dict[aa][0]
    stdev = aa_dict[aa][1]
    position = int(aa[0] - 1)
    aa_num = int(transform[aa[1]])
    try:
        fit_matrix[aa_num][position] = fitness
        dev_matrix[aa_num][position] = stdev
    except IndexError:
        pass
for i, col in enumerate(fit_matrix):
    for j, val in enumerate(col):
        if val == 0:
            fit_matrix[i][j] = 'NaN'

for i, col in enumerate(dev_matrix):
    for j, val in enumerate(col):
        if val == 0:
            dev_matrix[i][j] = 'NaN'




fit_matrix = np.ma.masked_invalid(fit_matrix)
dev_matrix = np.ma.masked_invalid(dev_matrix)

labels = np.arange(20)
plt.figure(1)
cmap = plt.get_cmap('YlGnBu')
cmap.set_under('white')
cmap.set_over('black')
cmap.set_bad(color = 'red', alpha = 0.5)
img = plt.imshow(dev_matrix, interpolation = 'nearest', cmap = cmap)
cb = plt.colorbar(img, extend = 'both' )       

plt.show()
outfile_name = (str(sys.argv[1][:-4]) + 'fit_matrix')


pic.dump(fit_matrix,open((str(outfile_name) + '.pkl'), 'wb'))

outfile_name = (outfile + '_dev_matrix')


pic.dump(dev_matrix,open((str(outfile_name) + '.pkl'), 'wb'))