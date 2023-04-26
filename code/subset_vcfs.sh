chr_num=$1
vcf=/vgipiper04/CCHC/TOPMed_postimpute_042022/chr${chr_num}.dose.vcf.gz
retained_samples=/data100t1/home/wanying/CCHC/lipidomics/input_docs/sample_list_for_bcftools.txt
output_fn=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/subset_chr${chr_num}.dose.vcf
echo "Processing chr${chr_num} with bcftools"
bcftools view -S ${retained_samples} -Ov ${vcf} > ${output_fn}

echo "bgzip file ${output_fn}"
bgzip ${output_fn}

