#!/bin/python
import os, sys
import argparse
from Bio import SeqIO
import cPickle as pickle
import seqmatch

'''Overall this scripts takes in the pkl file produced by the previous script for each time-point and day (TS) and checks the quality as 
well as the hamming distances of the barcodes. Moreover it counts how often you find each barcode (number of reads).
Parser: Provide the pkl produced by the fastq_index_parser.py The define an output name. The output file will again be a pkl file. Then 
provide the allele_pickle dictionary as reference. Finally define the hamming distance cutoff for the barcodes.
The script will then print out the stats for the direct, fuzzy and unmatched barcodes.'''

parser = argparse.ArgumentParser(description='Read in FASTQ file and indices, and saves pickle of timepoints with list of seqs with quality scores')
parser.add_argument('pkl_file', type=str,
                   help='fastq files to parse')
parser.add_argument('-o','--out_pickle', type=str, default="output_files/fuzzy_barcodes.pkl",
                   help='File where to save binary pickle')
parser.add_argument('--allele_pickle', type=str, help="allele pickle") 
parser.add_argument('--fuzzy_cutoff', type=int, default=None, help='Maximum hamming distance cutoff for matching barcodes [default: infinity]')

args = parser.parse_args()

#Listed the arguments in sys and renamed the 1st argument as the fastq file

output_dict = {}
barcodes_dict = pickle.load(open(args.allele_pickle, "rb"))
matched_index = pickle.load(open(args.pkl_file, "rb"))
final_counts = {}
barcode_keys = barcodes_dict.keys()

print "Pickles loaded."

for timepoint in matched_index.keys():
    tp_count_dict = { k:0 for k in barcode_keys  }
    total = 0
    direct = 0
    fuzzy = 0
    unmatched = 0
    print timepoint
    for barcode, quality in matched_index[timepoint]:
        total += 1

        if args.fuzzy_cutoff is None:
            (accepted_index, index_distance) = seqmatch.smart_match(barcode, barcode_keys)
        else:
            (accepted_index, index_distance) = seqmatch.smart_match(barcode, barcode_keys, max_dist=args.fuzzy_cutoff)
        
        if accepted_index is not None:
            if index_distance == 0:
                tp_count_dict[accepted_index] += 1
                direct += 1
            else:
                tp_count_dict[accepted_index] += 1
                fuzzy += 1
        else:
            unmatched += 1
    
    output_dict[timepoint] = tp_count_dict
    final_counts[timepoint] = (direct, fuzzy, unmatched, total)


for tp, fincounts in final_counts.items():
    print "Sample %s: %d exact matches, %d fuzzy matches, %d unmatched, %d total" % (tp, fincounts[0], fincounts[1], fincounts[2], fincounts[3])


pickle.dump(output_dict, open(args.out_pickle, 'w'))

print "done"

