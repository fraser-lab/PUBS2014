Module 0: Sub assembly
Script1:
01_sele_BC.py paired_end_read1.fasq > good_BC_reads.fastq

This script takes a raw fastq file (Illumina output) and checks each sequence to see if it matches the expected Ub-Library vector sequence after the N18 bar-code Input file should be the Read1 output file of a paired end Illumina read. The output is matched sequences in fastq format printed to the terminal. The script will also write a log file named "Script01_logfile.txt" by default

02_pair_reads.py good_BC_reads.fastq paired_end_read2.fastq 
output = pair_dict.pkl

This script takes the output fastq from 01_sele_BC.py and creates a dictionary keyed on the sequence sample ID. It then takes the full raw read2 fastq and associates the sample N18 bar-code with the Ub sequence from read2. The output is a dictionary keyed on the sample ID with values as a 2 item list. the first entry is the N18 bar-code (pair_dict[identity_key][0]). The second is a list of the Ub sequence from read2 (pair_dict[identity_key][1]). 

pair_dict.pkl
{‘SampleID’:[‘Barcode’, ‘Ub_sequence’], …}

03_sequences_assigned_to_barcode.py pair_dict.pkl 
output =  barcode_to_Ub.pkl

This script takes the output from 02_pair_reads.py and associates a given N18 barcode with all the related ubiquitin sequences. It then returns a dictionary that is keyed on the barcode with values of a list of all associated ubiquitin reads.

barcode_to_Ub.pkl
	{‘Barcode’: [‘Ub_sequence1’, ‘Ub_sequence2’, …] ...}
	
04_generate_consensus.py barcode_to_Ub.pkl
output = Allele_Dictionary.pkl

This script takes the output from 03_sequences_assigned_to_barcode.py and generates a consensus sequence from the list of Ub sequences associated with a given barcode. The mutant in the consensus sequence is identified and associated with the barcode. A barcode must be observed at least 3 times and the consensus sequence must contain only one mutation to be included. The output is a dictionary keyed on the barcode with a tuple value of (int(amino_acid_position), str(mutant_codon))

Allele_Dictionary.pkl
	{ “barcode” : (aa_position_in_Ub, Mutant_Codon)}
	“AGCTCTA” : (74, AUU)
“AGCCCTA” : (5, GCU)}
