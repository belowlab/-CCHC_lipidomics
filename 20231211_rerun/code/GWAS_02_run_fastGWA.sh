# Format of --pheno, --qcovar and --covar can be found here:
# https://gcta.freeforums.net/thread/510/limit-on-qcovars-fastgwa
# Only perform analysis in subst of samples specified by --keep
# Run this code in commandline as: ./GWAS_02_run_fastGWA.sh AC-OH class

lipid=$1
lipid_type=$2 # class or species
output_dir=/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_${lipid_type}
echo Processing ${lipid_type}

keep_list=/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/supporting_files/GWAS_files/samples_ids_for_fastGWA.train.txt

gcta64 --mbfile /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_imputed/plink_files.txt \
--grm-sparse /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_GRM/lipidomic_sp_grm \
--fastGWA-mlm \
--pheno /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/supporting_files/GWAS_files/lipid_${lipid_type}/${lipid}.pheno \
--keep ${keep_list} \
--thread-num 16 \
--out ${output_dir}/${lipid}