#!/bin/bash

# run_correlation_grid.sh
# Niels Hanson
# A script to do the coordinate the processing of SparCC results on Westgrid

# first calculate the correlation matrix on original dataset
input="EBPR_d_time.txt"
#input="fake_data.txt"

# downstream cor_to_network.py cutoffs
itr=20 # number of SparCC iterations
corr_thr=0.25 # SparCC min correlation threshold
num_boot=100 # Number of pval bootstap matrices

# network parameters
cor_cut_off=0.25 # min correlation cutoff
pval_cut_off=0.05 # max pvalue 

# westgrid directories
SparCC_dir="/home/nielsh/SparCC/"
qstat_dir="/home/nielsh/SparCC/qstat_dir"

# TODO clean up westgrid error output
# mkdir $qstat_dir

# python SparCC.py $input -i $itr --cor_file cor_out.txt -t $corr_thr
echo '#!/bin/bash' > temp.sh
echo '#PBS -l walltime=10:00:00' >> temp.sh
echo '#PBS -l pmem=2000mb' >> temp.sh
echo python ${SparCC_dir}SparCC.py ${SparCC_dir}${input} -i $itr --cor_file ${SparCC_dir}cor_out.txt -t ${corr_thr} >> temp.sh
# clean up
# echo rm ~/cov_mat_SparCC.out >> temp.sh
# echo mv ${SparCC_dir}temp.sh.e* ${qstat_dir} >> temp.sh
# echo mv ${SparCC_dir}temp.sh.o* ${qstat_dir} >> temp.sh
echo 1. Calcaulating initial correlation matrix...
echo Submitting: python SparCC.py $input -i $itr --cor_file cor_out.txt -t $corr_thr
echo ""
qsub temp.sh

# now shuffle the original data using makeBootstraps.py
echo 2. Generating suffled matrices...
echo Running: python ${SparCC_dir}MakeBootstraps.py ${SparCC_dir}${input} -n $num_boot -o "boot"
${SparCC_dir}MakeBootstraps.py ${SparCC_dir}${input} -n $num_boot -o "boot"
echo Done
echo ""

# now calculate the correlation on each shuffled "boot" matrix
echo 3. Calculating correlations on shuffled matrices...
for ((i=0; i<$num_boot; i++));
do
	echo '#!/bin/bash' > temp.sh
	echo '#PBS -l walltime=10:00:00' >> temp.sh
	echo '#PBS -l pmem=2000mb' >> temp.sh
	echo python ${SparCC_dir}SparCC.py ${SparCC_dir}boot_${i}.txt -i $itr --cor_file ${SparCC_dir}sim_cor_${i}.txt -t ${corr_thr}  >> temp.sh
	# clean up
	# echo rm ~/cov_mat_SparCC.out >> temp.sh
	# echo mv ${SparCC_dir}temp.sh.e* ${qstat_dir} >> temp.sh
	# echo mv ${SparCC_dir}temp.sh.o* ${qstat_dir} >> temp.sh
	
    echo Submitting: python SparCC.py boot_${i}.txt -i $itr --cor_file sim_cor_${i}.txt -t $corr_thr
	qsub temp.sh
done
echo Done
echo ""
exit 2

# downstream analyses
# echo 4. Calculating p-values of shuffled correlation matrices...
# echo python PseudoPvals.py cor_out.txt sim_cor $num_boot -o cor_out_pvals_one_sided.txt -t 'one_sided'
# python PseudoPvals.py cor_out.txt sim_cor $num_boot -o cor_out_pvals_one_sided.txt -t 'one_sided'
# echo Done
# echo ""
# 
# echo 5. Calculate network node and edge files...
# echo python cor_to_network.py -i cor_out.txt --pval cor_out_pvals_one_sided.txt --cor_cutoff $cor_cut_off --pval_cutoff $pval_cut_off
# python cor_to_network.py -i cor_out.txt --pval cor_out_pvals_one_sided.txt --cor_cutoff $cor_cut_off --pval_cutoff $pval_cut_off
# echo Done
# echo ""
# 
# # clean up intermediate files
# echo 6. Cleaning up extra files
# rm boot_*.txt
# rm sim_cor_*.txt
# echo Done
