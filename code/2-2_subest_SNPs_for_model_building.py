# This code extract SNPs from GWAS result using given p value threshold
# (Modified from 2_GWAS_create_plots.py)
# To run this code use:
# python 2_GWAS_create_plots.py --input input_file --output output_path
# 1. input path (GWAS result):
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS/*.fastGWA
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS/*.fastGWA
# 2. Output path:
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS_plots/
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_plots/

import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input') # input path and file name
parser.add_argument('-o', '--output') # output path (folder name)
args = parser.parse_args()

if not args.output.endswith('/'):
    args.output = args.output + '/'

# GWAS output format: CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P
lip = args.input.split('/')[-1].split('.fastGWA')[0]
print('\n\n################## Process:', lip, '##################')
print('# - Load GWAS result')
df = pd.read_csv(args.input, sep='\t')

threshold = 1e-3
print(f'# - Save varaints by pvalue (suggestive significant pval): pval<={threshold}')
df[df['P']<=threshold].to_csv(args.output+f'{lip}_SNPs_pval_{threshold}.txt', sep='\t', index=False)
print('################## DONE ##################')




