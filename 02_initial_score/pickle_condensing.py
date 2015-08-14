"""
Script that takes three pickles and condenses them into one pickle with a barcode as the key
and the value as a list of counts at the three time points for that day. 
the output will have the syntax of {Barcode1: [Time0Counts, Time1Counts, Time2Counts], Barcode2...}
"""
import cPickle as pic
import sys

T0 = sys.argv[1]
T1 = sys.argv[2]
T2 = sys.argv[3]
group = sys.argv[4]
day = sys.argv[5]

time0 = open(T0)
time1 = open(T1)
time2 = open(T2)

#Openning and loading the pickles for the three timepoints

pickle0 = pic.load(time0)
pickle1 = pic.load(time1)
pickle2 = pic.load(time2)

#cond will now be the condensed pickle
cond = {}

"""
To start the new dictionary, we grabbed all the keys with 15 or more counts from time 0 (this will help us avoid
some possible abnormalities in counts). This also adds the value associated with time 0
"""
for i in pickle0:
	cond[i] = [pickle0[i]]

#If the barcode in timepoint 1 matches with one of the barcodes from timepoint 0, append it to the value list

for i in pickle1:
	if i in cond.keys():
		cond[i].append(pickle1[i])	

#If the barcode didn't show up in timepoint 1, nothing was appended to the value. This checks and appends a 0. 

for key in cond.keys():
	if len(cond[key]) < 2:
		cond[key].append(0)	

#Same as before but with timepoint 2

for i in pickle2:
	if i in cond.keys():
		cond[i].append(pickle2[i])	

#Same as before but with timepoint 2

for key in cond.keys():
	if len(cond[key]) < 3:
		cond[key].append(0)			

#Save the file as a new pickle named whatever the args were that were present in the beginning (i.e. Dubstep 1)

pic.dump(cond, open("%s_%s.pkl" %(group, day), "wb"))