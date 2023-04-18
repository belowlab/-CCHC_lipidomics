chr_num=$1
tractor=/data100t1/home/wanying/downloaded_tools/Tractor/ExtractTracts.py
MSP_FILE=/vgipiper04/CCHC/local_ancestry/rfmix/CCHC_rfmix/output/CCHC_rfmix_chr${chr_num}
VCF_FILE=/vgipiper04/CCHC/TOPMed_postimpute_042022/chr${chr_num}.dose
output_dir=/data100t1/home/wanying/CCHC/lipidomics/output/snp_la_dosage/chr${chr_num}/

python ${tractor} \
--msp ${MSP_FILE} \
--vcf ${VCF_FILE} \
--output-path ${output_dir} \
--zipped \
--num-ancs 4
