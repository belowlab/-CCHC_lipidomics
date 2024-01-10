# Estimate heritability using filtered SNPs only.
# Execute the script in terminal

# example run: ./heritability_estimation_GREML_subset_SNPs.sh Sph-d18:1-

# Estimate heritability using filtered SNPs only.
# Execute the script in terminal

# example run: ./heritability_estimation_GREML_train_set_only.sh Sph-d18:1-

lipid=$1
lip_type='species'

output_dir=/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/heritability_by_GREML_species_subset_SNP_pval_1e-3
out=${output_dir}/${lipid}
grm_output_dir=${output_dir}/grm
pheno=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_${lip_type}_pheno/${lipid}.pheno
qcovar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-PC1-PC2-PC3-PC4-PC5.qcovar
covar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar
# grm=${grm_output_dir}/${lipid}
grm=/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/heritability_by_GREML_species_subset_SNP_pval_1e-3/grm/all_SNPs_combined_no_dup_no_multiallelic_species
snp_list=${output_dir}/${lipid}.snplist

# ########## Making a GRM from subset of SNPs in a family data set ##########
# Use *.snplist to only use a subset of SNPs
# gcta64 --mbfile /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/plink_files.txt \
# --keep /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/max_unrelated_set.indi.list \
# --extract ${snp_list} \
# --make-grm \
# --out ${grm} \
# --thread-num 16 \

# # ########## Creating an additional GRM from the GRM above (setting the off-diagonals that are < 0.05 to 0) ##########
# gcta64 --grm ${grm} --make-bK 0.05 --out ${grm}_bK --thread-num 16

# ########## Run a REML analysis on training set only, use test.indi.list to specify ##########
# This step ignore --extract flag
gcta64 --reml \
--keep /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/max_unrelated_set.indi.list \
--grm ${grm} \
--pheno ${pheno} \
--covar ${covar} \
--qcovar ${qcovar} \
--thread-num 16 \
--out ${out}