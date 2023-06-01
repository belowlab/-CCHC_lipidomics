# Use --keep to include training samples only
# Sample code
# ./run_fastGWA_lipd_species.sh AC-10:0-
lipid=$1
gcta64 --mbfile /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/plink_files.txt \
--grm-sparse /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/lipidomic_sp_grm \
--fastGWA-mlm \
--keep /data100t1/home/wanying/CCHC/lipidomics/input_docs/primus/primus_rel_3/list_for_fastGWA.txt \
--pheno /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_species_pheno/${lipid}.pheno \
--qcovar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-AGE2-PC1-PC2-PC3-PC4-PC5.qcovar \
--covar /data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar \
--save-fastGWA-mlm-residual \
--thread-num 8 \
--out /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/${lipid}

