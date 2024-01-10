# example run: ./heritability_estimation_GREML_train_set_only.sh Sph-d18:1-

lipid=$1
lip_type='species'

out=/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/heritability_by_GREML_${lip_type}/${lipid}
pheno=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_${lip_type}_pheno/${lipid}.pheno \
qcovar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-PC1-PC2-PC3-PC4-PC5.qcovar \
covar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar \
grm=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/lipidomic_grm 

# Running a REML analysis on training set only, use test.indi.list to specify
gcta64 --reml \
--keep /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/max_unrelated_set.indi.list \
--grm ${grm} \
--pheno ${pheno} \
--covar ${covar} \
--qcovar ${qcovar} \
--thread-num 32 \
--out ${out}