# Run a REML analysis on training set only, use a list to specify indlucded sample IDs
lipid=$1
lipid_type=$2 # class or species
covar=$3 # Whether to include covariates such as sex, age, PCs
output_dir=/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/heritability_estimation/${lipid_type}_no_covar
keep_list=/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/supporting_files/GWAS_files/samples_ids_for_fastGWA.train.txt
grm=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/lipidomic_grm

if [ ${covar} = 'false' ] || [ ${covar} = 'False' ] || [ ${covar} = 0 ]
then
  covars=''
else
  covars='--qcovar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-PC1-PC2-PC3-PC4-PC5.qcovar --covar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar'
fi

gcta64 --reml \
--keep ${keep_list} \
--grm ${grm} \
--pheno /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/supporting_files/GWAS_files/lipid_${lipid_type}/${lipid}.pheno \
--thread-num 32 \
--out ${output_dir}/${lipid} \
${covars}