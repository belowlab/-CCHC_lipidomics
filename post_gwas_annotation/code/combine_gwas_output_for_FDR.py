'''
Example usage:
python combine_gwas_output_for_FDR.py \
-i /data100t1/home/wanying/CCHC/lipidomics/output/lip_class_GWAS \
-o lipid_class_combined_bmi_age_age2_pc_sex.txt

python combine_gwas_output_for_FDR.py \
-i /data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS \
-o lipid_species_combined_bmi_age_age2_pc_sex.txt

'''

import pandas as pd
import numpy as np
from scipy import stats
import os
import datetime
import glob

import argparse

print(datetime.datetime.now().strftime('%Y-%m-%d'))
parser = argparse.ArgumentParser(description='Combine GWAS output for FDR correction')
parser.add_argument('-i', '--in_path')
parser.add_argument('-o', '--out_fn')

args = parser.parse_args()

path = args.in_path
output_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/lipdomics_GWAS_combined_FDR_corrected/'
# output_fn = 'lipid_class_combined_bmi_age_age2_pc_sex.txt'
output_fn = os.path.join(output_dir, args.out_fn)
lst_fn = glob.glob(path+'/*.fastGWA')
for i, fn in enumerate(lst_fn):
    # Cols in the output: CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P
    cols_to_use = ['SNP', 'P']
    df = pd.read_csv(fn, sep='\t')[cols_to_use]
    df['TRAIT'] = fn.split('.fastGWA')[0].split('/')[-1]
    print(f'\r# Files processed: {i+1}    ', end='', flush=True)
    if os.path.isfile(output_fn):
        df.to_csv(output_fn, mode='a', sep='\t', index=False, header=False)
    else:
        df.to_csv(output_fn, mode='a', sep='\t', index=False)
print('\n# DONE')