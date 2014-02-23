#!/bin/bash

# a script to do the coorelation analysis

# first calculate the correlation matrix on original dataset
input="EBPR_d_time.txt"
# input="fake_data.txt"
itr=20 # number of SparCC iterations
corr_thr=0.25 # SparCC min correlation threshold
num_boot=100 # Number of pval bootstap matrices

# downstream cor_to_network.py cutoffs
cor_cut_off=0.3 # min correlation cutoff
pval_cut_off=0.05 # max pvalue 

echo Calcaulating initial correlation matrix...
echo Running: python SparCC.py $input -i $itr --cor_file cor_out.txt -t $corr_thr
python SparCC.py $input -i $itr --cor_file cor_out.txt -t $corr_thr
echo Done

# now shuffle the original data using makeBootstraps.py
echo Generating suffled matrices...
echo Running: python MakeBootstraps.py $input -n $num_boot -o "boot"
python MakeBootstraps.py $input -n $num_boot -o "boot"
# python MakeBootstraps.py $input -n $num_boot -o "boot"
echo Done

# now calculate the correlation on each shuggled "boot" matrix
echo Calculating correlations on shuffled matrices...
for ((i=0; i<$num_boot; i++));
do
	echo Running: python SparCC.py boot_${i}.txt -i $itr --cor_file sim_cor_${i}.txt -t $corr_thr
	python SparCC.py boot_${i}.txt -i $itr --cor_file sim_cor_${i}.txt -t $corr_thr
done
echo Done

echo Calculating p-values from shuffed correaltion matrices...
echo python PseudoPvals.py cor_out.txt sim_cor $num_boot -o cor_out_pvals_one_sided.txt -t 'one_sided'
python PseudoPvals.py cor_out.txt sim_cor $num_boot -o cor_out_pvals_one_sided.txt -t 'one_sided'
echo Done

echo Calculate network node and edge files...
echo python cor_to_network.py -i cor_out.txt --pval cor_out_pvals_one_sided.txt --cor_cutoff $cor_cut_off --pval_cutoff $pval_cut_off
python cor_to_network.py -i cor_out.txt --pval cor_out_pvals_one_sided.txt --cor_cutoff $cor_cut_off --pval_cutoff $pval_cut_off
echo Done

# clean up intermediate files
echo Cleaning up extra files
rm boot_*.txt
rm sim_cor_*.txt
echo Done
