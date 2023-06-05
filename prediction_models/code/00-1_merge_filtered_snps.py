# Merge filtered SNPs of GWAS
# Remove duplications and sort by position, ref allele, alt allele
'''
Example code:
python 00-1_merge_filtered_snps.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_bmi_pval_1e-3/ \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_bmi_pval_1e-3/ \
--drop_multiallelic True

python 00-1_merge_filtered_snps.py \
--input_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/ \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/ \
--drop_multiallelic True
'''
import pandas as pd
import argparse
import os
import time
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/',
                    help='Input directory of filtered SNP files') # All files in this directory will be checked
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/') # Output file directory
parser.add_argument('--all_snps_fn',
                    default='all_SNPs_combined_no_dup_species.txt') # SNPs from all chromosomes, without duplication
parser.add_argument('--drop_multiallelic', default='False') # Whether to drop multiallelic sites
args = parser.parse_args()

if args.drop_multiallelic.upper()[0]=='T' or args.drop_multiallelic=='1':
    args.drop_multiallelic = True
else:
    args.drop_multiallelic = False

print('# Args used:')
print(args)

count = 0
cols = ['CHR', 'SNP', 'POS'] # Columns needed
for fn in os.listdir(args.input_dir):
    if count==0:
        df_all = pd.read_csv(args.input_dir+fn, sep='\t')[cols]
    else:
        df_tmp = pd.read_csv(args.input_dir+fn, sep='\t')[cols]
        df_all = df_all.merge(df_tmp, on=cols, how='outer')
    count += 1
    if count%250 == 0:
        print(f'{count} files processed')
    elif count%50==0:
        print('*', end='', flush=True)
    elif count%5 == 0:
        print('.', end='', flush=True)

print(f'\n# {count} files processed')
end_time = time.time()
print(f'# Merging finished in {(end_time - start_time)/60:.2f} min')

print('# Sort merged df')
df_all['REF'] = df_all['SNP'].apply(lambda x: x.split(':')[-2])
df_all['ALT'] = df_all['SNP'].apply(lambda x: x.split(':')[-1])
df_all.drop_duplicates(subset=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
df_all.sort_values(by=['CHR', 'POS', 'REF', 'ALT'], inplace=True)
if args.drop_multiallelic:
    df_all.drop_duplicates(subset=['CHR', 'POS'], keep=False, inplace=True)

print('\n\n# Merged df:', len(df_all))

df_all.to_csv(args.output_dir + args.all_snps_fn, sep='\t', index=False)
print(f'# Saving finished in {(time.time() - end_time)/60:.2f} min')
print('# DONE')
