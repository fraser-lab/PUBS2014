#!/bin/python
import os, sys
import argparse
from Bio import SeqIO
import cPickle as pickle
import seqmatch

# Constants
READ_CONSTANT = "AGATCGG"

'''Overall this scripts extracts the indices from the fastq files and save all the barcodes and quality scores according to time-points and days (TS) into 
a pkl file (dictionaries). 
Parser: Define here which fastq file you want to read in (day and lane). Moreover you need to provide a txt file with all indice barcode sequences
 (here called  indices.txt. Then define the output name. This will be a pkl file which then can be read as input into the next script pkl_barcode_parser.py
 In addition you can provide the hamming distance for the indices (--index_cutoff) and the hamming distance for the constant region (--const_cutoff).
 Finally this scripts prints the stats for direct match and fuzzy match of indices as well as the fuzzy matched constant regions.'''
parser = argparse.ArgumentParser(description='Read in FASTQ file and indices, and saves pickle of timepoints with list of seqs with quality scores')
parser.add_argument('fastq_files', type=str, nargs='+',
                   help='fastq files to parse')
parser.add_argument('--indices', type=str, default=None,
                   help='tab-delimited file with two columns (timepoint first column, index seq second column)')
parser.add_argument('-o','--out_pickle', type=str, default="barcodes_dict.pkl",
                   help='File where to save binary pickle')
parser.add_argument('--write_fastq', type=str, default=None, help='write filtered sequences to FASTQ file')
parser.add_argument('--index_cutoff', type=int, default=2, help='distance cutoff to use for index')
parser.add_argument('--const_cutoff', type=int, default=2, help='distance cutoff to use for constant region')

args = parser.parse_args()

if aesrgs.indic is None or not os.path.isfile(args.indices):
    print "Please specify an indices file with --indices"
    parser.print_usage()
    sys.exit()

#Listed the arguments in sys and renamed the 1st argument as the fastq file
indices_file = args.indices
fastq_files = args.fastq_files
barcodes_dict = {}
indices_dict = {}
indices_list = []
index_dists = {}
constant_dists = {}

with open(indices_file,"rU") as indices_handle:
    for line in indices_handle:
        timepoint,seq_index = line.strip().split("\t")
        indices_dict[seq_index] = timepoint
indices_list = indices_dict.keys()

if args.write_fastq is not None:
    output_handle = open(args.write_fastq, "w")

for fastq_file in fastq_files:
    file_handle = open(fastq_file, "rU")

    #create a parser to parse file handle with seqio
    parser = SeqIO.parse(file_handle,"fastq")

    for record in parser:
        quality = record.letter_annotations["phred_quality"]
        badscores = [i for i in quality if i<20]
        if len(badscores) < 3:
            record_index = record.description.strip()[-6:]
            (accepted_index, index_distance) = seqmatch.smart_match(record_index,indices_list,max_dist=args.index_cutoff)
            if accepted_index is not None:
                seq = record.seq
                constant = str(seq[-7:])
                (accepted_constant, constant_distance) = seqmatch.smart_match(constant, [READ_CONSTANT,], max_dist=args.const_cutoff)
                if accepted_constant is not None:
                    time_point = indices_dict[accepted_index]
                    
                    if time_point not in constant_dists:
                        constant_dists[time_point] = {}
                    if constant_distance not in constant_dists[time_point]:
                        constant_dists[time_point][constant_distance]=0
                    constant_dists[time_point][constant_distance]+=1

                    if time_point not in index_dists:
                        index_dists[time_point] = {}
                    if index_distance not in index_dists[time_point]:
                        index_dists[time_point][index_distance]=0
                    index_dists[time_point][index_distance]+=1

                    quality = record.letter_annotations["phred_quality"][:18]
                    barcode = str(seq[:18].reverse_complement())
                    
                    if time_point not in barcodes_dict:
                        barcodes_dict[time_point] = []
                    barcodes_dict[time_point].append( (barcode,quality) )
                    
                    if args.write_fastq is not None:
                        record.description = "%s%s" % (record.description[:-7],accepted_index)
                        record.seq = "%s%s" % (seq[:18], accepted_constant)
                        SeqIO.write(record,output_handle, "fastq")

    if args.write_fastq is not None:
        output_handle.close()

for tp, bc_list in barcodes_dict.items():
    index_fuzzy_count = 0
    constant_fuzzy_count = 0
    
    if tp in index_dists:
        for dist,dist_count in [(x,y) for x,y in index_dists[tp].items() if x > 0]:
            index_fuzzy_count += dist_count

    if tp in constant_dists:
        for dist,dist_count in [(x,y) for x,y in constant_dists[tp].items() if x > 0]:
            constant_fuzzy_count += dist_count


    print "Sample %s: %d barcodes (%d fuzzy index, %d fuzzy constant)" % (tp, len(bc_list), index_fuzzy_count, constant_fuzzy_count)


pickle.dump(barcodes_dict, open(args.out_pickle, 'w'))

