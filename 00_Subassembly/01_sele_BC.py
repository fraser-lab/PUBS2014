'''
-This script takes a raw fastq file (Illumina output) and checks each sequence to see if it matches the expected Ub-Library vector sequence after the N18 bar-code
-Input file should be the Read1 output file of a paired end Illumina read
-The output is matched sequences in fastq format printed to the terminal. Please > this output to a file of your choice. The script will also write a log file named "Script01_logfile.txt"
by default. Feel free to change this name at line 28
Cheers
'''

import Bio
from Bio import SeqIO
import sys

in_file_1 = open(sys.argv[1], 'rU') #Raw sequencing reads (.fastq)

# the sequence we are looking for 22 bp downstream of the bar-code
constant = 'GGCAGGACTCAGGGCCGCTGCA' 
count_good = 0.0 # how many reads match the constant region
count_total = 0.0 # total number of reads analyzed

for record in SeqIO.parse(in_file_1, 'fastq'):
	count_total += 1
	if constant ==  str(record.seq[18:40]): # Are the 22 bp downstream of the N18 bar-code the constant region?
		count_good	+= 1
		print record.format('fastq') # print to a new .fastq with only the good reads
in_file_1.close()	

log_file =  open('Script01_logfile.txt', 'w') # writes a logfile with total number of reads, the number of reads picked as 'good' and the percentage of 'good reads'

log_file.write('Good reads = %d \n' % (count_good) )
log_file.write('Total reads =  %d \n' % (count_total))
log_file.write('Percent good reads =  %f \n' % (count_good/count_total))
log_file.close()