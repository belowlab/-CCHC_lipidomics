lipid=$1
gcta64 --mbfile /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/plink_files.txt \
	--grm-sparse /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/lipidomic_sp_grm \
	--fastGWA-mlm \
	--pheno /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_class_pheno/${lipid}.pheno \
	--qcovar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-PC1-PC2-PC3-PC4-PC5.qcovar \
	--covar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar \
	--save-fastGWA-mlm-residual \
	--thread-num 16 \
	--out /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS_noadj_BMI_AGE2/${lipid}
