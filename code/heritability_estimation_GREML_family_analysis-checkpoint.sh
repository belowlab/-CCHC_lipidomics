lipid=$1

# GREML in family data
out=/data100t1/home/wanying/CCHC/lipidomics/output/heritability_by_GREML_species/${lipid}
pheno=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/lipid_species_pheno/${lipid}.pheno \
qcovar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/AGE_AT_VISIT-PC1-PC2-PC3-PC4-PC5.qcovar \
covar=/data100t1/home/wanying/CCHC/lipidomics/input_docs/pheno_covar_files/GENDER.covar \
# Multiple GRMs
mgrm_file=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/mgrm.txt

# Running a REML analysis with two GRMs
gcta64 --reml \
--mgrm ${mgrm_file} \
--pheno ${pheno} \
--covar ${covar} \
--qcovar ${qcovar} \
--thread-num 32 \
--out ${out}
