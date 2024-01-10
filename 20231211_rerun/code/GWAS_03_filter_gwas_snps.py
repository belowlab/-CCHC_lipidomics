# This code extract SNPs from fastGWAS result using given p value threshold and MAF
# (Modified from /data100t1/home/wanying/CCHC/lipidomics/code/2_GWAS_create_plots.py)
# To run this code use:
# python GWAS_03_filter_gwas_snps.py --input input_file --output output_path --pval 1e-3 --maf 0.05
# - If filter by MAF is not needed
# python GWAS_03_filter_gwas_snps.py --input input_file --output output_path --pval 1e-3

# 1. input path (fastGWA results):
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_class/*.fastGWA
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_species/*.fastGWA
# 2. Output path:
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS_plots/
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_plots/

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input') # input path and file name
parser.add_argument('-o', '--output_dir') # output path (folder name)
parser.add_argument('-t', '--pval', type=float, default=1e-3) # Threshold of pvalue to filter SNPs
parser.add_argument('-m', '--maf', type=float, default=-1) # Threshold of MAF, -1 means no filtering
args = parser.parse_args()

# GWAS output format: CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P
lip = args.input.split('/')[-1].split('.fastGWA')[0]
print('\n\n#', '#'*20, 'Process:', lip, '#'*20)
print('# - Load GWAS result')
df = pd.read_csv(args.input, sep='\t')

print(f'# - Filter variants by pvalue (suggestive significant pval): pval<={args.pval}')
if not os.path.isdir(os.path.expanduser(args.output_dir)): # Create the directory if not exist
    print('# - Output path does not exist, creating one:', os.path.expanduser(args.output_dir))
    os.mkdir(os.path.expanduser(args.output_dir))
    
if args.maf==-1:
    mask = df['P']<=args.pval
    out_fn = os.path.join(args.output_dir, f'{lip}_SNPs.pval_{args.pval}.txt')
else:
    args.maf = float(args.maf)
    print(f'# - Filter variants by MAF>={args.maf}')
    out_fn = os.path.join(args.output_dir, f'{lip}_SNPs.pval_{args.pval}.maf_{args.maf}.txt')
    mask = (df['P'] <= args.pval) & (df['AF1']>=args.maf) & (df['AF1']<=1-args.maf)
    
df[mask].to_csv(out_fn, sep='\t', index=False)
print('#', '#'*20, 'DONE', '#'*20)
