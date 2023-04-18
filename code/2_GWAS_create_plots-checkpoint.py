# To run this code use:
# python 2_GWAS_create_plots.py --input input_file --output output_path
# 1. input path (GWAS result):
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS/*.fastGWA
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS/*.fastGWA
# 2. Output path:
#  (1) Lipid class: /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS_plots/
#  (2) Lipid species: /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_plots/

import pandas as pd
import matplotlib as plt
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/')
from manhattan_plot import manhattan_plot
from QQplot_v6 import qqplot
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input') # input path and file name
parser.add_argument('-o', '--output') # output path
parser.add_argument('-ht', '--heritability') # Heritability
args = parser.parse_args()

if not args.output.endswith('/'):
    args.output = args.output + '/'

# GWAS output format: CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P
lip = args.input.split('/')[-1].split('.fastGWA')[0]
print('\n\n################## Process:', lip, '##################')
print('# - Load GWAS result')
df = pd.read_csv(args.input, sep='\t')
print('# - Create Manhattan plot')
fig_mh, ax_mh = manhattan_plot(data=df, pval='P', position='POS', chromosome='CHR',
                               sig_pval=0.05/(10**6), title=f'{lip}, heritability={args.heritability}')
fig_mh.savefig(f'{args.output}{lip}_{args.heritability}_manhattan_plot.jpeg')

print('# - Create QQ plot')
fig_qq, ax_qq, lambda_original, lambda_novel = qqplot(df, output=f'{args.output}{lip}_{args.heritability}_QQ_plot.jpeg',
                                                      p_value_column_title='P',
                                                      title=f'{lip}, heritability={args.heritability}')

print('# - Save varaints by pvalue (suggestive significant pval): pval<=1e-5')
df[df['P']<=1e-5].to_csv(args.output+f'{lip}_suggestive_sig_SNPs.txt', sep='\t', index=False)
print('################## DONE ##################')




