"""
-The input is a dictionary (pickle object) keyed on the id with values as a 2 item list. the first entry is the N18 bar-code. The second is a list of the Ub sequence from read2.  
-The output is a dictionary (pickle) keyed on the barcode with a list of sequences as the value
"""

import sys
import cPickle as pic

pair_dict = pic.load( open( sys.argv[1], "rb" ) ) #input dictionary from command line, this is the output of 02

barcode_to_sequence_dict = {} 

for key in pair_dict:
	barcode = pair_dict[key][0] 
	ub_read_sequence = pair_dict[key][1] 
	try: #checks if first part of value list (barcode) has already been seen
		barcode_to_sequence_dict[barcode].append(ub_read_sequence)
	except KeyError:
		barcode_to_sequence_dict[barcode] = [ub_read_sequence]

pic.dump(barcode_to_sequence_dict, open('barcode_to_Ub.pkl', 'wb'))