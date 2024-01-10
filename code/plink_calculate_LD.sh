# Calcualte LD of a given set of SNPs
# Input bfile has to be a single chromosme for now
# Example input
# input: /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/chr1
# snp list: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS/all_SNPs_combined_no_dup_pval_1e-8_chr1.txt

chr=$1

input=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/chr${chr}
snp_list=/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS/all_SNPs_combined_no_dup_pval_1e-8_chr${chr}.txt

echo "input files: ${input}, ${snp_list}"

# Need to set ld-window to a large number to avoid filtering of SNPs
# Reference: https://zzz.bwh.harvard.edu/plink/ld.shtml#ld2
plink --bfile ${input} \
--r2 \
--ld-window-kb 2000 \
-ld-window 99999 \
--ld-window-r2 0 \
--ld-snp-list ${snp_list} \
--out ./plink/chr${chr}