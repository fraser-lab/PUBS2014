import subprocess
import sys

position = sys.argv[1]
aa =  sys.argv[2]


files = ['05_heatmaps/DMSO_ave_to_Barcode_FitScore_Caffeine_outliers_removed_matrix.pkl', '05_heatmaps/DMSO_ave_to_Barcode_FitScore_DTT_outliers_removed_matrix.pkl', '05_heatmaps/DMSO_ave_to_Barcode_FitScore_HU_outliers_removed_matrix.pkl']

for i in files:
	to_call = 'python output_score.py %s %s %s' %(i, position, aa)
	print to_call
	subprocess.call(to_call, shell = True)