import pandas as pd
import os
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/')
from extract_snps_for_large_vcf_v2_xopen import find_variants
# from multiprocessing import Pool

'''# Try result of the simple model: trait ~ sex + age + snp + PC1-5 + grm
result_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_noadj_BMI_AGE2_snps_pval_1e-5/'
input_dir = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/'
output_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/prediction_models/elastic_net/'

# Load SNPs with pval<1e-5
for fn in os.listdir(result_dir):
    if fn.endswith('.txt'):
        lipid = fn.split('_suggestive_sig_SNPs.txt')[0]
        print(f'# Processing lipid {lipid}')
        df_snps = pd.read_csv(result_dir+fn, sep='\t')
        df_snps.sort_values(by=['CHR', 'POS'], inplace=True)

        for chr_num, df in df_snps.groupby('CHR'):
            print(f'# - Processing chr{chr_num}')
            input_fn = f'max_unrelated_set_chr{chr_num}.vcf.gz'
            output_fn = f'training_{lipid}_chr{chr_num}.vcf'
            with Pool(30) as p:
                p.starmap(find_variants, [(df_snps['POS'], output_dir+output_fn, input_dir+input_fn, 'POS', 8)])'''


# Use unrealted samples in training and testing
# Save others for validation
input_dir = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/'
output_dir = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/'
lip_type = 'species'

# Load SNPs to be extracted
print('# Load SNPs to be extracted')
all_snps_fn = '/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3/all_SNPs_combined_no_dup_species.txt'
df_snps = pd.read_csv(all_snps_fn, sep='\t')

for chr_num, df in df_snps.groupby('CHR'):
    input_fn = f'max_unrelated_set_chr{chr_num}.vcf.gz'
    output_fn = f'{lip_type}_chr{chr_num}.vcf'
    print(f'\n# Process chr {chr_num}')
    print(f'# - Subset VCF saved to {output_fn}')
    # with Pool(30) as p:
    #     p.starmap(find_variants, [(df_snps['POS'], output_dir+output_fn, input_dir+input_fn, 'POS', 8)])
    find_variants(lst_pos = df['POS'],
                  output_fn = output_dir+output_fn,
                  input_fn = input_dir+input_fn,
                  input_col_name = 'POS',
                  threads = 8)
