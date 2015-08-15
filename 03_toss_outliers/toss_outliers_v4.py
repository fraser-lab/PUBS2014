import cPickle as pic 
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Input: dictionary of fitness scores corresponding to a codon or mutation
	# Either:
	# {barcode:[fitness, stde, pval, rval]}
	# {(pos,codon):[fitness, stde, pval, rval]}
bc_mut_dict = pic.load(open('pkl/allele_dic_with_WT.pkl',"rb"))
translate = pic.load(open('pkl/translate.pkl',"rb"))
# Run
# cd to directory
# python toss_outliers_v3.py pickle.pkl sort_type cutoff
# where sort_type should be either 'barcode' or 'codon'
# where cutoff should be an int, probably 4 for barcode and 0 for codon

# CREATES A MATRIX OF DICTIONARIES, VALUE=FITNESS, SIZE AND KEYS VARY BY SORT_TYPE
def create_fitness_mat(filename, sort_type='codon'):

	# Calling all pkls
	A2N = pic.load(open('pkl/aa_to_num_sorted.pkl',"rb"))
	bc_mut_dict = pic.load(open('pkl/allele_dic_with_WT.pkl',"rb"))
	translate = pic.load(open('pkl/translate.pkl',"rb"))
	fitness_dict = pic.load(open(filename, 'rb'))
	codon_ypos = pic.load(open('pkl/codon_ypos.pkl', 'rb'))

	# Save WT Ub aa sequence, obsolete
	# wtUb = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG*'
	# wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTTGA'
	
	# Establish number of aa positions
	p = 77
	wt_fit = []

	# If fitness was determined by barcode
	if sort_type == 'barcode':
		q = 64

		# Make Npos x Ncodon matrix of dictionaries, {barcode:fitness}
		fitness_mat = [[{} for i in range(q)] for j in range(p)]
		for barcode in fitness_dict.keys():
			pos, codon = bc_mut_dict[barcode]
			if not codon == 'WT':
				rna = codon.replace('T','U')
				aa = translate[rna]
				# wtaa = wtUb[pos-1]
				# if not aa == wtaa:
				if not np.isnan(fitness_dict[barcode][0]):
					fitness_mat[pos-1][codon_ypos[codon]][barcode] = (fitness_dict[barcode][0], fitness_dict[barcode][-1])
				else:
					print fitness_dict[barcode],
			else:
				wt_fit.append(fitness_dict[barcode][0])


	# If fitness was determined by codon
	elif sort_type == 'codon':
		q = 21
		# Make Npos x Naa matrix of dictionaries, {(pos,codon):fitness}
		fitness_mat = [[{} for i in range(q)] for j in range(p)]
		for barcode in fitness_dict.keys():
			pos, codon = bc_mut_dict[barcode]
			if not codon == 'WT':
				rna = codon.replace('T','U')
				aa = translate[rna]
				# wtaa = wtUb[pos-1]
				# if not aa == wtaa:
				if not np.isnan(fitness_dict[barcode][0]):
					fitness_mat[pos-1][A2N[aa]][(pos,codon,barcode)] = (fitness_dict[barcode][0], fitness_dict[barcode][-1])
				else:
					print pos, codon,fitness_dict[(pos,codon)],
			else:
				wt_fit.append(fitness_dict[barcode][0])

	wt_std = np.std(wt_fit)
	wt_mean = np.mean(wt_fit)
	for fit in wt_fit:
		if abs(fit-wt_mean) > 2*wt_std:
			wt_fit.remove(fit)

	wt_mean = np.mean(wt_fit)

		
	return fitness_mat, wt_mean

# COMPUTES A HEATMAP OF THE DIFFERENCE IN VARIANCE 
def var_heatmap(fitness_mat, sort_type, wt_mean, cutoff=4):

	A2N = pic.load(open("pkl/aa_to_num_sorted.pkl", "rb"))
	codon_ypos = pic.load(open("pkl/codon_ypos.pkl", "rb"))
	translate = pic.load(open("pkl/translate.pkl","rb"))
	codon_dic = pic.load(open('pkl/codon_ypos.pkl', "rb"))
	pathname = os.path.split(filename)

	# Determine shape of matrix & make empty lists
	p, q = np.shape(fitness_mat)

	before_mat = [[[] for i in range(q)] for j in range(p)]
	after_mat = [[[] for i in range(q)] for j in range(p)]
	clean_list = []
	dirty_list = []

	# Iterate through matrix, get outliers for each
	# separate into clean and dirty dictionaries for heatmap
	# append to lists for ultimate pkl output
	for pos in range(p):
		for c in range(q):
			dict_barcodes = fitness_mat[pos][c]
			if len(dict_barcodes) > cutoff:
				clean, dirty = get_outliers(dict_barcodes, wt_mean)
			else:
				clean = dict_barcodes
				dirty = {}

			for barcode in clean.keys():
				# if not np.isnan(fitness):
				clean_list.append(barcode)
				after_mat[pos][c].append(clean[barcode])
			for fitness in dict_barcodes.values():
				# if not np.isnan(potentialvalue):
				before_mat[pos][c].append(fitness)
			for barcode in dirty.keys():
				# if not np.isnan(potentialvalue):
				dirty_list.append(barcode)

	#COMPUTE VARIANCES BY POS AND AA/CODON
	fit_var_before = np.zeros((p,q))
	fit_var_after = np.zeros((p,q))
	for pos in range(p):
		for c in range(q):
			if len(after_mat[pos][c]) > 1:
				fit_var_after[pos,c]=np.std(np.array(after_mat[pos][c]))
			if len(before_mat[pos][c]) > 1:
				fit_var_before[pos,c]=np.std(np.array(before_mat[pos][c]))	

	# ESTABLISH MIN/MAX TO NORMALIZE HEATMAPS
	minv = min([np.amin(fit_var_before)] + [np.amin(fit_var_after)])
	maxv = max([np.amax(fit_var_before)] + [np.amax(fit_var_after)])
	# w = [x/10.0 for x in range(-1,4)]

	# IMPORTANT FOR ACCURATE TICKS ON HEATMAP
	n = [.5]+list(np.arange(4.5,76.5,5))
	m = [1]+list(np.arange(5,77,5))

	fig = plt.figure()
	fig.subplots_adjust(hspace=.3)

	# PLOT HEATMAPS
	if sort_type == 'barcode':
		p = 77
		q = 21

		# Plot heatmaps
		ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
		ax.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax.tick_params(axis='both', which='minor', labelsize=8)
		ax.tick_params(axis='both', which='major', labelsize=5)
		ax.set_xticks(n, minor=True)
		ax.set_xticklabels(m, minor=True)
		ax.set_yticks(np.array(codon_ypos.values())+.5)
		ax.set_yticklabels(codon_ypos.keys())
		ax.set_xlim([0,p])
		ax.set_title('Before tossing outliers', fontsize=10)
		ax.set_aspect('equal')
		plt.pcolormesh((fit_var_before.T), cmap='Spectral', vmin=minv, vmax=maxv)
		cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
		cbar.ax.tick_params(labelsize=10) 

		ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
		ax2.tick_params(axis='both', which='minor', labelsize=8)
		ax2.tick_params(axis='both', which='major', labelsize=5)
		ax2.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax2.set_xticks(n, minor=True)
		ax2.set_xticklabels(m, minor=True)
		ax2.set_yticks(np.array(codon_ypos.values())+.5)
		ax2.set_yticklabels(codon_ypos.keys())
		ax2.set_xlim([0,p])
		ax2.set_title('After tossing outliers', fontsize=10)
		ax2.set_aspect('equal')
		plt.pcolormesh((fit_var_after.T), cmap='Spectral', vmin=minv, vmax=maxv)
		cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
		cbar2.ax.tick_params(labelsize=10) 

		ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
		ax3.tick_params(axis='both', which='minor', labelsize=8)
		ax3.tick_params(axis='both', which='major', labelsize=5)
		ax3.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax3.set_xticks(n, minor=True)
		ax3.set_xticklabels(m, minor=True)
		ax3.set_yticks(np.array(codon_ypos.values())+.5)
		ax3.set_yticklabels(codon_ypos.keys())
		ax3.set_xlim([0,p])
		ax3.set_title('Subtraction of before-after', fontsize=10)
		ax3.set_aspect('equal')
		plt.pcolormesh((fit_var_before.T-fit_var_after.T), cmap='Spectral')
		cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
		cbar3.ax.tick_params(labelsize=10) 
		
		heatmaps_d = os.path.dirname(filename)+'/heatmaps'
		if not os.path.exists(heatmaps_d):
			os.makedirs(heatmaps_d)

		plt.savefig(pathname[0]+'/heatmaps/'+pathname[1]+'_fitness_variance_by_codon.png')

	if sort_type == 'codon':
		p = 77
		q = 21

		fig = plt.figure()
		fig.subplots_adjust(hspace=.3)

		ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
		ax.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax.tick_params(axis='both', which='minor', labelsize=8)
		ax.tick_params(axis='both', which='major', labelsize=5)
		ax.set_xticks(n, minor=True)
		ax.set_xticklabels(m, minor=True)
		ax.set_yticks(np.array(A2N.values())+.5)
		ax.set_yticklabels(A2N.keys())
		ax.set_xlim([0,p])
		ax.set_title('Before tossing outliers', fontsize=10)
		ax.set_aspect('equal')
		plt.pcolormesh((fit_var_before.T), cmap='Spectral', vmin=minv, vmax=maxv)
		cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
		cbar.ax.tick_params(labelsize=10) 

		ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
		ax2.tick_params(axis='both', which='minor', labelsize=8)
		ax2.tick_params(axis='both', which='major', labelsize=5)
		ax2.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax2.set_xticks(n, minor=True)
		ax2.set_xticklabels(m, minor=True)
		ax2.set_yticks(np.array(A2N.values())+.5)
		ax2.set_yticklabels(A2N.keys())
		ax2.set_xlim([0,p])
		ax2.set_title('After tossing outliers', fontsize=10)
		ax2.set_aspect('equal')
		plt.pcolormesh((fit_var_after.T), cmap='Spectral', vmin=minv, vmax=maxv)
		cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
		cbar2.ax.tick_params(labelsize=10) 

		ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
		ax3.tick_params(axis='both', which='minor', labelsize=8)
		ax3.tick_params(axis='both', which='major', labelsize=5)
		ax3.set_xticklabels('')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ax3.set_xticks(n, minor=True)
		ax3.set_xticklabels(m, minor=True)
		ax3.set_yticks(np.array(A2N.values())+.5)
		ax3.set_yticklabels(A2N.keys())
		ax3.set_xlim([0,p])
		ax3.set_title('Subtraction of before-after', fontsize=10)
		ax3.set_aspect('equal')
		plt.pcolormesh((fit_var_before.T-fit_var_after.T), cmap='Spectral')
		cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
		cbar3.ax.tick_params(labelsize=10)

		heatmaps_d = os.path.dirname(filename)+'/heatmaps'
		if not os.path.exists(heatmaps_d):
			os.makedirs(heatmaps_d)

		plt.savefig(pathname[0]+'/heatmaps/'+pathname[1]+'_fitness_variance_by_aa.png')

	return before_mat, after_mat, clean_list, dirty_list

def get_outliers(fitness_dict, wt_mean, thresh=1.5):
			fit_tup = [(k, fitness_dict[k][0], fitness_dict[k][1]) for k in fitness_dict]
			scores = [k[1]-wt_mean for k in fit_tup]
			#calculate MAD
			# med = np.median(scores)
			# mad = 1.4826*np.median(abs(np.subtract(scores, med)))

			# numerator = abs(scores)- float(med)
			# mad_scores = numerator/mad

			clean_dict = {}
			outliers = {}
			med = np.median(scores)
			q75, q25 = np.percentile(scores,[75, 25])
			IQR = q75-q25


			for i in range(len(fit_tup)):
				key = fit_tup[i][0]
				codon = key[1]
				rna = codon.replace('T','U')
				fitness = fit_tup[i][1]
				std_err = fit_tup[i][2]
				lower = q25 - thresh*IQR
				upper = q75 + thresh*IQR
			
				if fitness >= lower and fitness <= upper:
					clean_dict[key] = fitness
				else:
					outliers[key] = fitness

			return clean_dict, outliers

if __name__ == "__main__":
	filename = sys.argv[1]
	sort_type = sys.argv[2]
	cutoff = int(sys.argv[3])

	# make matrix of dictionaries, with fitness as values
	fitness_mat, wt_mean = create_fitness_mat(filename, sort_type)

	# compute and output heatmap by aa
	before, after, clean_list, dirty_list = var_heatmap(fitness_mat, sort_type, wt_mean, cutoff)

	# PRINT NUMBER OF CLEAN AND DIRTY BARCODES FOUND, UNCOMMENT BELOW TO ALSO
	# PRINT NUMBER OF FITNESS SCORES THROWN OUT
	# fitness_dict = pic.load(open(filename, 'rb'))
	# x = len(fitness_dict) - len(clean_list) - len(dirty_list)
	print len(clean_list), '\tclean barcodes'
	print len(dirty_list), '\tdirty barcodes'
	# print x, '\tWTs thrown out'

	# OUTPUT DICT
	pathname = os.path.split(filename)
	output_dict = {'clean_barcodes':clean_list, 'dirty_barcodes':dirty_list}
	out_d = os.path.dirname(filename)+'/output_new/'
	if not os.path.exists(out_d):
		os.makedirs(out_d)
	pic.dump(output_dict, open(out_d + pathname[1] + '_clean_list.pkl', 'w'), protocol=2)

	with open(out_d+'list_counts.txt', 'a') as tossfile:
		# txt file: pkl_file clean_counts dirty_counts
		newline = [pathname[1], str(len(clean_list)), str(len(dirty_list))]
		catline = '\t'.join(newline)
		tossfile.write(catline+'\n')







