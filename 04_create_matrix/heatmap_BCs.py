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

matrix = np.zeros((21,77))
table = string.maketrans('T', 'U')
big_list= [[[] for i in range(77)] for i in range(21)]
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

#codons
# for key in dictionary.keys():
#     aa_num = key[0]
#     codon = key[1]
#     codon = string.translate(codon, table)
#     aa_val = int(TRANSLATE_num[codon])
    # cords = (aa_val-1, int(aa_num)-1)  
    # print cords
    # matrix[cords[0]][cords[1]] = dictionary[key]

#aap
codon_dict = {}
for key in dictionary:
    codon = allele_dict[key]
    fit_score = dictionary[key][0]
    if codon in codon_dict.keys():
        codon_dict[codon].append(fit_score)
    else:
        codon_dict[codon] = [fit_score]
for codon in codon_dict.keys():
    codon_dict[codon] =  average_list(codon_dict[codon])


for key in codon_dict.keys():
    if key[1] != 'WT':
        aa_num = key[0]
        codon = key[1]
        codon = string.translate(codon, table)
        aa = str(TRANSLATE[codon])
        aa_val = int(transform[aa])
        cords = (aa_val, (int(aa_num)-1))  
        big_list[cords[0]][cords[1]].append(codon_dict[key])


min_tracker = []

for i in range(len(big_list)):
    for j in range(len(big_list[0])):
        try:
            min_tracker.append(min(big_list[i][j]))
        except ValueError:
            pass
        try:
            matrix[i][j] = average_list(big_list[i][j])
        except ZeroDivisionError:
            matrix[i][j] = 'NaN'


matrix = np.ma.masked_invalid(matrix)
labels = np.arange(20)
print labels
plt.figure(1)
cmap = plt.get_cmap('YlGnBu')
cmap.set_under('white')
cmap.set_over('black')
cmap.set_bad(color = 'red', alpha = 0.5)
img = plt.imshow(matrix, interpolation = 'nearest', vmax = 0.2, vmin = -1.0, cmap = cmap)
cb = plt.colorbar(img, extend = 'both' )       
#img.set_yticks(np.array(range(len(labels))) + 0.5)
#img.set_yticklabels(labels) 
#plt.show()
outfile_name = (str(sys.argv[1][:-4]) + '_matrix')

plt.savefig((str(outfile_name) + '.png'))

pic.dump(matrix,open((str(outfile_name) + '.pkl'), 'wb'))

# # count_of_barcodes = defaultdict(int)

# # for pair in dictionary:
# #   count_of_barcodes[dictionary[pair][0]] +=1


# # sorted_x = sorted(count_of_barcodes.iteritems(), key=operator.itemgetter(1))

# # for x in sorted_x:
# #   if x[0][0] == "T":
# #       print x[0],x[1]
# # print sorted_x


# # for bc in count_of_barcodes:
# #   print bc, count_of_barcodes[bc]






# # pylab.hist(coun