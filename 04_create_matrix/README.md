Module 4: Create matrix
heatmap_BCs.py Barcode_scores_outliers_removed.pkl Allele_Dictionary.pkl
output = Barcode_scores_outliers_removed_fit_matrix.pkl & Barcode_scores_outliers_removed_dev_matrix.pkl

This script takes the Barcode scores and averages them to the amino acid level and calculates the standard devation of the barcodes mapping to a given amino acid substitution. It then outputs these scores as a heatmap and the scores and SDs as numpy matrix pkls. 

Barcode_scores_outliers_removed_fit_matrix.pkl
masked_array(data = [21X76 matrix containing fitness scores for each aa substitution])

[[-- -0.631791406846642 -0.5397724613430753 ..., -0.3530569856099873
 -0.4209070611721436 --]
 [-- -0.04155544432657808 -- ..., -- -0.4672829444990562
 -0.015863306980341812]
 [-- -0.06881685222913404 -0.3000996826283508 ..., -0.09038257104760622
 -0.5060198247122988 --]
 ...,
 [-- -0.2136845391374962 -0.5846954623554699 ..., -0.5201027046981986 -- --]
 [-- -0.037103840372513595 -0.6445621511224743 ..., -0.6398418859301449
 -0.5009308293341608 --]
 [-- -0.03813077561871194 -0.5946959237696324 ..., -0.5684471534238328
 -0.3959407495759722 --]]

Barcode_scores_outliers_removed_dev_matrix.pkl
masked_array(data = [21X76 matrix containing the standard deviation of the barcode fitness scores for each aa substitution])

[[-- -- 0.07810572453203178 ..., 0.03738608565402174 -- --]
 [-- -- -- ..., -- -- --]
 [-- 0.09548287851233793 -- ..., -- -- --]
 ...,
 [-- -- -- ..., 0.02225009725522345 -- --]
 [-- 0.04200403236126699 -- ..., -- -- --]
 [-- 0.05500820984278576 -- ..., 0.10087705678557354 0.11792320714288289 --]]