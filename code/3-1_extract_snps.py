import pandas as pd
import os
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/')
from extract_snps_for_large_vcf_v2_xopen import find_variants
from multiprocessing import Pool

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', default=None) # Input file directory
parser.add_argument('--output_dir', default=None) # Output file directory
parser.add_argument('--result_dir', default=None) # Result file directory
parser.add_argument('--gwas_summary', default=None) # a single gwas summary file file to be processed, will scan the entire input directory if not provided
parser.add_argument('--type', default='training') # Training or validation

args = parser.parse_args()

# Try result of the simple model: trait ~ sex + age + snp + PC1-5 + grm
if not args.result_dir:
    args.result_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_noadj_BMI_AGE2_snps_pval_1e-5/'
if not args.input_dir:
    args.input_dir = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/validation/'
if not args.output_dir:
    args.output_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/prediction_models/elastic_net/validation/'
if not args.gwas_summary:
    args.gwas_summary = os.gwas_summary(result_dir)
else:
    args.gwas_summary = [args.gwas_summary]

# Load SNPs with pval<1e-5
for fn in args.gwas_summary:
    if fn.endswith('.txt'):
        lipid = fn.split('_suggestive_sig_SNPs.txt')[0]
        print(f'# Processing lipid {lipid}')
        df_snps = pd.read_csv(args.result_dir+fn, sep='\t')
        df_snps.sort_values(by=['CHR', 'POS'], inplace=True)
        
        for chr_num, df in df_snps.groupby('CHR'):
            print(f'# - Processing chr{chr_num}')
            if args.type=='training':
                input_fn = f'max_unrelated_set_chr{chr_num}.vcf.gz'
            else:
                input_fn = f'{args.type}_set_chr{chr_num}.vcf.gz'
            output_fn = f'{args.type}_{lipid}_chr{chr_num}.vcf'
            with Pool(30) as p:
                p.starmap(find_variants, [(df_snps['POS'], args.output_dir+output_fn, args.input_dir+input_fn, 'POS', 8)])


