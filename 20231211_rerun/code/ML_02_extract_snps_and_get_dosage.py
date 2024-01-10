# Train regression model (elastic net or ridge) using filtered GWAS SNPs
# Modified from ~/CCHC/lipidomics/code/3-1_extract_snps_and_get_dosage.py
# First create a master file contains filtered SNPs from all lipid traits

'''
Call script as:
python ML_02_extract_snps_and_get_dosage.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs \
--input_fn max_unrelated_set_chr*.vcf.gz  \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train \
--output_fn lipid_species_chr*.pval_0.001_maf_0.05.vcf \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_class_filter_by_pval_1e-3_MAF_5e-2 \
--all_snps_fn all_SNPs.pval_0.001_maf_0.05.txt
--chr_range 1-22

python ML_02_extract_snps_and_get_dosage.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs \
--input_fn max_unrelated_set_chr*.vcf.gz  \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_class/tmp \
--output_fn lipid_class_chr*.pval_0.001_maf_0.05.vcf \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_class_filter_by_pval_1e-3 \
--all_snps_fn all_SNPs.pval_0.001.txt \
--chr_range 3

# For test set
python ML_02_extract_snps_and_get_dosage.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_test_vcfs \
--input_fn test_set_chr*.vcf.gz  \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test \
--output_fn lipid_species_chr*.pval_0.001_maf_0.05.vcf \
--result_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_class_filter_by_pval_1e-3_MAF_5e-2 \
--all_snps_fn all_SNPs.pval_0.001_maf_0.05.txt
'''

# ################# TODO #################
# Need to change code to account for multiallelic sites. Cannot just look by SNP positions

import pandas as pd
import sys

sys.path.append('/data100t1/home/wanying/lab_code/utils/')
import os
from extract_snps_for_large_vcf_tabix import find_variants
from get_dosage_from_vcf_compressed import get_dosage
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs',
                    help='Input vcf file directory')
parser.add_argument('--input_fn',
                    default='max_unrelated_set_chr*.vcf.gz ')  # Input file format, use * to represent chromosome number
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train')  # Output file directory
parser.add_argument('--output_fn',
                    default='subset_chr*.vcf')  # output file name of subsetted vcf, use * to represent chromosome number
parser.add_argument('--result_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/')  # Filtered SNP result file directory
parser.add_argument('--all_snps_fn',
                    default='all_SNPs_combined_no_dup_no_multiallelic_species.txt')  # SNPs from all chromosomes
parser.add_argument('--chr_range', default='1-22') # Range of chromosome or a signle chr number to process. Default to check chr1 to chr22
args = parser.parse_args()

# Sanity checks
print('# Sanity checks')
print('# - Check input VCF files:')
all_input_file_exist = True  # Flag

# Start and end chromosome number to be processed
chr_start, chr_end = int(args.chr_range.split('-')[0]), int(args.chr_range.split('-')[-1])
for i in range(chr_start, chr_end+1):  # Need to check from chromosome 1 to 22 (or chr_start to chr_end)
    if not os.path.isfile(os.path.join(os.path.expanduser(args.input_dir),
                                       args.input_fn.replace('*', str(i)))):
        print(f"#\t - chr{i}: {args.input_fn.replace('*', str(i))} does not exist")
        all_input_file_exist = False
    else:
        print(f'#\t - CHR{i}: PASS')
if not all_input_file_exist:
    print('\n# DONE')
    exit()

print('# - Check merged and filtered SNP file: ', end='')
if not os.path.isfile(os.path.join(os.path.expanduser(args.result_dir), args.all_snps_fn)):
    print(f'\n# - Input file does not exist. Exit')
    print('\n# DONE')
    exit()
else:
    print('PASS')

print('# - Check output directory: ', end='')
if not os.path.isdir(os.path.expanduser(args.output_dir)):
    print(f'\n# - Output directory does not exist, create one at: {os.path.expanduser(args.output_dir)}')
    os.mkdir(os.path.expanduser(args.output_dir))
else:
    print('Output directory exists. PASS')

# Load SNPs to be extracted
print('\n# Load SNPs to be extracted')
df_snps = pd.read_csv(os.path.join(args.result_dir, args.all_snps_fn), sep='\t')
df_snps.sort_values(by=['CHR', 'POS', 'A1', 'A2'], inplace=True)

for chr_num, df in df_snps.groupby('CHR'):
    if chr_num>=chr_start and chr_num<=chr_end:
        input_fn = args.input_fn.replace('*', str(chr_num))
        output_fn = args.output_fn.replace('*', str(chr_num))
        print(f'\n# Process chr {chr_num}: {input_fn}')
        print(f'# - Subset VCF saved to {os.path.join(args.output_dir, output_fn)}')

        # For this version use chr1:1234 to locate variants,
        # but might need to change and match chromosome format in vcf, such as 1:1234 or CHR1:1234, etc
        lst_pos = [f'chr{x[0]}:{x[1]}' for x in list(df[['CHR', 'POS']].itertuples(index=False, name=None))]
        find_variants(lst_pos=lst_pos,
                      output_fn=os.path.join(args.output_dir, output_fn),
                      input_fn=os.path.join(args.input_dir, input_fn),
                      keep_multiallelic=False, bgzip=True)
        # Then extract dosage
        try:
            print(f'# - Get dosage')
            output_dosage_fn = os.path.join(args.output_dir, output_fn) + '.dosage'
            get_dosage(vcf_fn=os.path.join(args.output_dir, output_fn) + '.gz',
                       output_fn=output_dosage_fn)
            print('# - bgzip and tabix index dosage file')
            cmd = f'bgzip {output_dosage_fn}; tabix -b 2 -e 2 -f {output_dosage_fn}.gz'
            subprocess.run(cmd, shell=True)
        except:
            print('# - Get dosage failed, do it manually:', os.path.join(args.output_dir, output_fn) + '.gz')
print('\n# DONE')
