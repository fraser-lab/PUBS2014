#!/bin/python
import os, sys
import argparse
from Bio import SeqIO
import cPickle as pickle
import seqmatch

'''This script finally converts the previous dictionaries into the requested order. One dictionary for each TS with barcodes sequence as key
and the number of reads as values. 
Parser: As input provide the pkl files (there should be three for Oct16-lane1, Oct16-lane2 and Oct17-lane3). Then define an output directory and 
find a basename. That means each pkl for each TS will be named as follows: basename_TSxx where xx runs from 1 to 30. Again, you need to 
provide the allele_pickle dictionary as reference.'''

parser = argparse.ArgumentParser(description='Read in FASTQ file and indices, and saves pickle of timepoints with list of seqs with quality scores')
parser.add_argument('pkl_file', type=str, nargs='+',
                   help='counts pkl files to parse')
parser.add_argument('-o','--out_dir', type=str, default="output_files",
                   help='Output directory for saving pickles [default: outdir]')
parser.add_argument('-b','--pkl_basename', type=str, default="final_counts",
                   help='basename for saving final output pickles [default: final_counts]')
parser.add_argument('--allele_pickle', type=str, help="allele pickle")
args = parser.parse_args()

#Listed the arguments in sys and renamed the 1st argument as the fastq file
allele_dict = pickle.load(open(args.allele_pickle, "rb"))
barcode_keys = allele_dict.keys()

barcode_counts_dicts = []
for pkl_file in args.pkl_file:
    barcode_counts_dicts.append(pickle.load(open(pkl_file, "rb")))

all_indices = set([])
for barcode_counts_dict in barcode_counts_dicts:
    indices = barcode_counts_dict.keys()
    print indices
    all_indices.update(indices)
    
for sample_index in sorted(all_indices):
    output_pkl = "%s/%s_%s.pkl" % (args.out_dir, args.pkl_basename, sample_index)
    final_counts_dict = {k:0 for k in barcode_keys}
    for barcodes_counts_dict in barcode_counts_dicts:
        if sample_index in barcodes_counts_dict:
            sample_barcodes_counts = barcodes_counts_dict[sample_index]
            for barcode, count in barcodes_counts_dict[sample_index].items():
                final_counts_dict[barcode] += count
    pickle.dump(final_counts_dict,open(output_pkl,'w'))
    
