import pandas as pd
import os
import datetime
print(datetime.datetime.now().strftime('%Y-%m-%d'))


# #################### Change here ####################
output_fn = 'SNP2gene_mapped.txt'
window = 0 # Eg. 500kb +/- gene region
# Load reference and SNP file
df_gene_ref = pd.read_csv('../supporting_files/ensembl_ref_GRCh38.103.txt', sep='\t')
df_combined = pd.read_csv('../input/combined_snps_1e-8_lipid_species_gwas_noadj_bmi_age2.txt', sep='\t')
# #################### Do not change below this line ####################

# Map SNP to gene region
c = 0 # Count
not_mapped_gene = 0 # Number of SNPs without mapped gene
mapped_genes = []
total_num_snps = len(df_combined)
for i in range(total_num_snps):
    # Get chromosome number and SNP position
    chr_num, pos = df_combined.iloc[i]['CHR'], df_combined.iloc[i]['POS']
    # Look for gene in the reference file
    mask = (df_gene_ref['chr']==f'{chr_num}') & ((df_gene_ref['start']-window<=pos) & (df_gene_ref['end']+window>=pos))
    df_result = df_gene_ref[mask]
    if len(df_result) > 0:
        mapped_genes.append(','.join(df_result['gene']))
    else:
        mapped_genes.append(None) # No gene mapped
        not_mapped_gene += 1
    c += 1
    print(f'\r# Number of SNPs processed: {c}/{total_num_snps}    ', end='', flush=True)
print(f'\n# Number of SNPs without mapped gene: {not_mapped_gene}')
print('# DONE')

df_combined['mapped_gene'] = mapped_genes
df_combined.to_csv(output_fn, sep='\t', index=False)