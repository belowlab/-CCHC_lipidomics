from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
import argparse
import warnings
import datetime
print('Last run:', datetime.datetime.now().strftime('%Y-%m-%d'))
warnings.filterwarnings(action='ignore')

# ---------------------- Help functions ----------------------
def get_doasge(dosage_fn, lst_snps):
    '''
    Param:
     - dosage_fn: name of dosage file to be check against (single chromosome only)
     - lst_snps: a list of SNP positions to be searched for (single chromosome only)
    Return:
     - sample_ids: smaple IDs
     - dosage_matrix: dosage of given SNPs as a numpy array. Fill with NA if a SNP is not found
    '''
    
    with open(dosage_fn) as fh:
        line = fh.readline().strip() # Take sample IDs from header line
        tmp = line.split()
        indx_dosage = tmp.index('FORMAT') + 1 # Get index of sample IDs and dosage values
        indx_pos = tmp.index('POS') # Index of SNP position
        sample_ids = tmp[indx_dosage:] # Genotype IDs
        dosage_matrix = [] # Store dosage values in a numpy matrix. Lines are SNPs, columns are individuals
        
        line = fh.readline().strip()
        snp_pos = lst_snps.pop(0) # Check from the first element
        count = 0
        print('\t', end='')
        while line != '':
            # Scan through dosage file to get dosage of GWAS snps
            tmp = line.split()
            cur_pos = tmp[indx_pos]
            
            if float(cur_pos) == float(snp_pos): # Find a match
                dosage = tmp[indx_dosage:]
                dosage_matrix += dosage
                if len(lst_snps) > 0:
                    snp_pos = lst_snps.pop(0)
                else:
                    break
                line = fh.readline().strip()
                count += 1
            elif float(cur_pos) > float(snp_pos):
                # print(dosage_fn, cur_pos) # For testing !!!!
                dosage = [np.nan] * len(sample_ids) # SNP not found in dosage file, fill dosage with NAs
                dosage_matrix += dosage
                if len(lst_snps) > 0:
                    # If current position in dosage file is already larger than SNP pos
                    # Does not need to read in the next line
                    snp_pos = lst_snps.pop(0) # Check next SNP
                else:
                    # Does not need to continue reading dosage file when the SNP list is empty
                    break
            else:
                # Keep reading in the next line if SNP pos is smaller than current pos
                line = fh.readline().strip()
                count += 1
            
            if count%1000000==0:
                print(f'{count} lines processed', flush=True)
                print('\t', end='')
            elif count%20000==0:
                print('.', end='', flush=True)
    print(f'{count} lines processed')            
    return sample_ids, np.array(dosage_matrix).reshape(-1, len(sample_ids))

# Load dosage of all SNPs with p val<10-3 from GWAS
def load_all_dosage(gwas_snp_fn: str,
                    gwas_snp_dir: str='',
                    dosage_dir: str='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train',
                    dosage_fn: str='species_chr*.vcf.gz.dosage'):
    '''
    Get doage of all SNPs (GWAS pval<1e-3) from single-chrosmosome dosage files of a given lipid
    Params:
        - gwas_snp_dir: directory to GWAS SNPs
        - gwas_snp_fn: file name of GWAS SNPs
        - dosage_dir: Subsetted dosage file: species_chr*.vcf.gz.dosage
        - dosage_fn: file name of subset dosage files (by chromosome).
                    Replace chromosome number with '*', such as 'species_chr*.vcf.gz.dosage'
    Return:
        - df_gwas_snp: a dataframe of GWAS SNPs
        - dosage_all: A numpy array of doage. Each row is a SNP, each column is a subject
    '''
    # Check if file exists
    if gwas_snp_dir.endswith('/'): gwas_snp_dir = gwas_snp_dir[:-1] # Remove last slash
    if not os.path.isfile(f'{gwas_snp_dir}/{gwas_snp_fn}'):
        print(f'# ERROR: GWAS SNP file not find: {gwas_snp_dir}/{gwas_snp_fn}\n# END')
        exit()
        
    lip_name = gwas_snp_fn.split('_')[0]
    print('# Processing lipid:', lip_name)

    # print(f'# Load GWAS SNPs for current lipid')
    df_gwas_snp = pd.read_csv(f'{gwas_snp_dir}/{gwas_snp_fn}', sep='\t').sort_values(by=['CHR', 'POS'])
    # print(f'# - Number of SNPs loaded: {len(df_gwas_snp)}')

    print('\n# Get dosage of GWAS SNPs to include in regression models')
    print('# - Checking by chromosome:')

    dosage_all = '' # A numpy array to store dosage from all chromosome
    start_time = datetime.datetime.now() # Time execution time
    for chr_num, df in df_gwas_snp.groupby(by='CHR'):
        # dosage_fn = f'species_chr{chr_num}.vcf.gz.dosage'
        print(f'#  chr{chr_num}')
        sample_ids, dosage_matrix = get_doasge(f"{dosage_dir}/{dosage_fn.replace('*', str(chr_num))}", list(df['POS']))
        # lst_df_dosage.append(pd.DataFrame(data=dosage_matrix, columns=sample_ids, index=df['POS']))
        if len(dosage_all) == 0: # if dosage array is empty
            dosage_all = dosage_matrix
        else:
            dosage_all = np.append(dosage_all, dosage_matrix, axis=0)
        # break
    end_time = datetime.datetime.now()
    print(f'# - Checking finished in {(end_time-start_time).total_seconds()}s')
    print('-' * 50)
    return df_gwas_snp, dosage_all.astype('float64')

# ---------------------- End of help functions ----------------------


# ################# Process args #################
parser = argparse.ArgumentParser(description='Fit elastic net regression with 10 fold cross-validation',
                                 epilog='Text at the bottom of help')

# ################# Load lipidomic data #################
print('# Load lipidomic data (lipid species)')
fn_lipid = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_measures/lipid_species.txt'
df_lipid = pd.read_csv(fn_lipid, sep='\t')
print(f"# - data loaded from {fn_lipid.split('/')[-1]}: shape {df_lipid.shape}")

# Re-order lipidomic data so that sample IDs match the order in genotype file
fn_id_mapping = '/data100t1/home/wanying/CCHC/doc/samples_IDs/202211_merged_RNA_lipid_protein_genotype_mapping_and_availability.txt'
df_id_mapping = pd.read_csv(fn_id_mapping,
                            sep='\t').dropna(subset=['genotype_ID',
                                                     'lipidomic']).drop_duplicates(subset='genotype_ID')[['LABID', 'genotype_ID']]

print(f'\n# Load genotype IDs for matching (only need to read the first line of dosage file)')
dosage_dir = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train'
fn_genotype = f'{dosage_dir}/species_chr22.vcf.gz.dosage'
with open(fn_genotype) as fh:
    df_genotype_id = pd.DataFrame(fh.readline().strip().split()[9:], columns=['genotype_ID'])

print(f'# - Organize sample IDs so that their orders match in lipidomics data and dosage file')
df_lipid = df_genotype_id.merge(df_id_mapping.merge(df_lipid.drop_duplicates(subset='Sample ID'),
                                                    left_on='LABID',
                                                    right_on='Sample ID'), on='genotype_ID')
print(f'# - Final processed lipidomic data: {len(df_lipid)}')


# ################# Load GWAS snps of each lipid and run regression #################
# dosage_all: each row contains doages of a single SNP across all individuals
# !! Lip species PI(15-MHDA_20:4)\PI(17:0_20:4) is missing
gwas_snp_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3' # GWAS SNPs with p value<1e-3
output_file = f"{datetime.datetime.now().strftime('%Y%m%d-%M:%S')}_lip_species_elasticnet_params.txt" # Save coefficients, alpha and l1 ratios of selected model for each lipid
output_fh = open(output_file, 'w')
output_fh.write('lipid\talpha\tl1_ratio\tcoefficients\n') # write header line

count = 0 
for lip in df_lipid.columns[4:]:
    gwas_snp_fn = f"{lip.replace('(', '-').replace(')', '-').replace(' ', '_').replace('/', '-')}_SNPs_pval_0.001.txt"
    if os.path.isfile(f'{gwas_snp_dir}/{gwas_snp_fn}'):
        lip_name = gwas_snp_fn.split('_')[0] # Modified lipid name
        # Get SNPs and dosage
        print(f'\n# Load GWAS SNPs for current lipid: {lip_name}')
        df_gwas_snp,dosage_all = load_all_dosage(gwas_snp_dir = gwas_snp_dir,
                                                 gwas_snp_fn = gwas_snp_fn,
                                                 dosage_dir = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train',
                                                 dosage_fn = 'species_chr*.vcf.gz.dosage')
        print(f'# - Number of SNPs loaded: {len(df_gwas_snp)}')
        
        print('# Run Elastic net regression')
        # lipid level, INVed
        y = inverse_normal_transformation(df_lipid[lip])
        # print(y.shape)

        start_time = time.time()
        # regr = ElasticNet(alpha=0.5, max_iter=10000, random_state=0)
        # regr = ElasticNet(alpha=0.5, random_state=0)
        # alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        alphas = [0.25, 0.5, 0.75, 1]
        # l1_ratio = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1] # It is recommanded to put more values close to 1 (i.e. Lasso) and less close to 0 (i.e. Ridge)
        # regr = ElasticNetCV(cv=10, random_state=0, n_jobs=32, alphas=alphas, l1_ratio=l1_ratio)
        regr = ElasticNetCV(cv=10, random_state=0, n_jobs=32, alphas=alphas) # Try 11 ratio=0.5 first
        regr.fit(dosage_all.T, y)

        end_time = time.time()
        print(f'# - Model fitting finised in {(end_time - start_time):.4f}s')
        output_fh.write(f"{lip}\t{regr.alpha_}\t{regr.l1_ratio_}\t{','.join(str(x) for x in regr.coef_)}\n")
        # alpha\tl1_ratio\tcoefficients\n'
        # break
    else:
        print(f'# - Warning: {lip} not found')
    count += 1
    print(f'# #################### {count} lipid processed ####################')
output_fh.close()
