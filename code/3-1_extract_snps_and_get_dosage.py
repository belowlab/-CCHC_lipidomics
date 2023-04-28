'''
Call script as:
python 3-1_extract_snps.py --input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/validation/ \
--input_fn validation_set_chr*.vcf.gz \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/validation/ \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3/ \
--all_snps_fn all_SNPs_combined_no_dup_species.txt \
--lip_type species
'''
import pandas as pd
import os
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/')
from extract_snps_for_large_vcf_v2_xopen import find_variants
from get_dosage_from_vcf_compressed import get_dosage
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/validation/') # Input file directory
parser.add_argument('--input_fn',
                    default='validation_set_chr*.vcf.gz') # Input file format, use * to represent chromosom number
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/validation/') # Output file directory
parser.add_argument('--result_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3/') # Result file directory
parser.add_argument('--all_snps_fn',
                    default='all_SNPs_combined_no_dup_species.txt') # SNPs from all chromosomes, without duplication
parser.add_argument('--lip_type', default='species', choices=['class', 'species']) # Lipid class or species
args = parser.parse_args()

# Load SNPs with pval<1e-3 of model: trait ~ sex + age + age2 + BMI + snp + PC1-5 + grm
# Use unrealted samples in training and testing
# Save others for validation
# input_dir = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/'
# output_dir = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/'
# lip_type = 'species'

# Load SNPs to be extracted
print('# Load SNPs to be extracted')
df_snps = pd.read_csv(f'{args.result_dir}{args.all_snps_fn}', sep='\t')

for chr_num, df in df_snps.groupby('CHR'):
    # input_fn = f'max_unrelated_set_chr{chr_num}.vcf.gz'
    input_fn = args.input_fn.replace('*', str(chr_num))
    output_fn = f'{args.lip_type}_chr{chr_num}.vcf'
    print(f'\n# Process chr {chr_num}: {input_fn}')
    print(f'# - Subset VCF saved to {output_fn}')
    # with Pool(30) as p:
    #     p.starmap(find_variants, [(df_snps['POS'], output_dir+output_fn, input_dir+input_fn, 'POS', 8)])
    find_variants(lst_pos = df['POS'],
                  output_fn = args.output_dir+output_fn,
                  input_fn = args.input_dir+input_fn,
                  input_col_name = 'POS',
                  threads = 8)
    # Then extract dosage
    print(f'# - Get dosage')
    get_dosage(vcf_fn=args.output_dir+output_fn)
