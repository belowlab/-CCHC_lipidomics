'''
Call script as:
python 3-1_extract_snps_and_get_dosage.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/ \
--input_fn max_unrelated_set_chr*.vcf.gz  \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/ \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/ \
--all_snps_fn all_SNPs_combined_no_dup_no_multiallelic_species.txt \
--lip_type species


# For test set
python 3-1_extract_snps_and_get_dosage.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/test/ \
--input_fn test_set_chr*.vcf.gz \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/test/ \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/ \
--all_snps_fn all_SNPs_combined_no_dup_no_multiallelic_species.txt \
--lip_type species

'''

# ################# TODO #################
# Need to change code to account for multiallelic sites. Cannot just look by SNP positions

import pandas as pd
import os
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/')
from extract_snps_for_large_vcf_v2_xopen_multiallelic import find_variants_multiallelic
from get_dosage_from_vcf_compressed import get_dosage
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/training_max_unrelated_sampels_3rd_degree/',
                    help='Input vcf file directory')
parser.add_argument('--input_fn',
                    default='max_unrelated_set_chr*.vcf.gz') # Input file format, use * to represent chromosome number
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/') # Output file directory
parser.add_argument('--result_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/') # Filtered result file directory
parser.add_argument('--all_snps_fn',
                    default='all_SNPs_combined_no_dup_no_multiallelic_species.txt') # SNPs from all chromosomes, without duplication
parser.add_argument('--lip_type', default='species', choices=['class', 'species']) # Lipid class or species
args = parser.parse_args()

if not args.input_dir.endswith('/'): args.input_dir = args.input_dir+'/'
if not args.output_dir.endswith('/'): args.output_dir = args.output_dir+'/'

# Load SNPs with pval<1e-3 of model: trait ~ sex + age + snp + PC1-5 + grm
# Use unrelated samples in training and testing
# Save others for validation

# Load SNPs to be extracted
print('# Load SNPs to be extracted')
df_snps = pd.read_csv(f'{args.result_dir}{args.all_snps_fn}', sep='\t')
df_snps.sort_values(by=['CHR', 'POS', 'REF', 'ALT'], inplace=True)

for chr_num, df in df_snps.groupby('CHR'):
    input_fn = args.input_fn.replace('*', str(chr_num))
    output_fn = f'{args.lip_type}_chr{chr_num}.vcf'
    print(f'\n# Process chr {chr_num}: {input_fn}')
    print(f'# - Subset VCF saved to {args.output_dir}{output_fn}')
    find_variants_multiallelic(lst_pos_ref_alt = list(df[['POS', 'REF', 'ALT']].itertuples(index=False, name=None)),
                               output_fn = args.output_dir+output_fn,
                               input_fn = args.input_dir+input_fn,
                               input_col_names = ['POS', 'REF', 'ALT'],
                               threads = 8, verbose=True)
    # Then extract dosage
    print(f'# - Get dosage')
    try:
        get_dosage(vcf_fn=args.output_dir+output_fn)
    except:
        print('# - Get dosage failed, do it manually:', args.output_dir+output_fn)
