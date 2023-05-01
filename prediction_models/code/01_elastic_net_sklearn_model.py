from sklearn.linear_model import ElasticNetCV
import pandas as pd
import numpy as np
import os
import time
import argparse
import sys
# sys.path.append('/data100t1/home/wanying/lab_code/utils')
# from rank_based_inverse_normal_transformation import inverse_normal_transformation
import warnings
import datetime
warnings.filterwarnings(action='ignore')

'''
Example call:
# python 01_elastic_net_sklearn_model.py --output lipid_species_l1_0.5_500-599.txt --output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params --lipid_range 500 --range_window 100 --n_alphas 10

python 01_elastic_net_sklearn_model.py --output lipid_species_l1_0.5_100-199.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params \
--lipid_range 100 \
--range_window 100 \
--n_alphas 100

# (Per AlexP) To run in base environment and avoid multi-threading conflict
# Run script as: OMP_NUM_THREADS=1 python my_script.py, for example:

OMP_NUM_THREADS=1 python 01_elastic_net_sklearn_model.py --output lipid_species_l1_0.5_5-104_100_alpha_CV.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/100alphas \
--lipid_range 5 \
--range_window 100 \
--n_alphas 100

'''
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
        - start_time: start time of loading dosage
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
    start_time = time.time() # Time execution time
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
    end_time = time.time()
    print(f'# - Checking finished in {(end_time-start_time):.4f}s')
    print('-' * 50)
    return start_time, df_gwas_snp, dosage_all.astype('float64')

# ---------------------- End of help functions ----------------------


# ################# Process args #################
parser = argparse.ArgumentParser(description='Fit elastic net regression with 10 fold cross-validation')
parser.add_argument('-o', '--output', type=str,
                           help='Output file to  save alpha, l1_ratio and coefficients of chosen model')
parser.add_argument('--output_dir', type=str, help='Output directory. Defualt is current directory', default='.')
parser.add_argument('--lipid_range', type=int, default=-1,
                    help='Define a subset of lipids to run. Default -1 ie. run all lipids')
parser.add_argument('--range_window', type=int, default=100,
                    help='Define a window of lipids to run. Default is 100, will run from lipid_range to lipid_range+range_window-1')
parser.add_argument('--n_alphas', type=int, default=100,
                    help='Define how many alphas to test in CV. Dafault is 10. JTI used 100 as defined in R glmnet()')

args = parser.parse_args()
args.output = f"{args.output}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"
if args.output_dir.endswith('/'): args.output_dir = args.output_dir[:-1]
print('# Run starts:', datetime.datetime.now().strftime('%Y-%m-%d'))
print('# Output file is', f'{args.output_dir}/{args.output}')
print(f'# Cross validation on {args.n_alphas} alphas')
if args.lipid_range==-1:
    print('# Run all lipids. Ignore range_window')
else:
    print(f'# Run lipids from index {args.lipid_range} to {args.lipid_range + args.range_window - 1}')

# ################# Load lipidomic data #################
print('# Load lipidomic data (lipid species)')
# fn_lipid = '/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_measures/lipid_species.txt'
fn_lipid = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted'
df_lipid = pd.read_csv(fn_lipid, sep='\t')
print(f"# - data loaded from {fn_lipid.split('/')[-1]}: shape {df_lipid.shape}")

# Re-order lipidomic data so that sample IDs match the order in genotype file
fn_id_mapping = '/data100t1/home/wanying/CCHC/doc/samples_IDs/202211_merged_RNA_lipid_protein_genotype_mapping_and_availability.txt'
df_id_mapping = pd.read_csv(fn_id_mapping,
                            sep='\t').dropna(subset=['genotype_ID',
                                                     'lipidomic']).drop_duplicates(subset='lipidomic')[['LABID', 'genotype_ID']]

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
# Save coefficients, alpha and l1 ratios of selected model for each lipid
output_fh = open(f'{args.output_dir}/{args.output}', 'w')
output_fh.write('lipid\talpha\tl1_ratio\tbest_r2\tcoefficients\n') # write header line
output_fh_lip_pred = open(f'{args.output_dir}/{args.output}.pred', 'w') # Save predicted values of each lipid using best fitted model
output_fh_lip_pred.write('Lipid'+'\t'+'\t'.join([val for val in df_lipid['Sample ID']])+'\n') # write header line

count = 0
# get list of lipids to be fitted
if args.lipid_range != -1:
    if args.lipid_range>=len(df_lipid.columns[4:]):
        print(f'# ERROR: lipid range {args.lipid_range} is out of range')
        exit()
    else:
        try:
            lst_lipids = df_lipid.columns[3+args.lipid_range : 4+args.lipid_range+args.range_window]
        except:
            # If window is too large, just run to the end of lipid list
            lst_lipids = df_lipid.columns[3 + args.lipid_range:]
else:
    lst_lipids = df_lipid.columns[3:]

for lip in lst_lipids:
    gwas_snp_fn = f"{lip.replace('(', '-').replace(')', '-').replace(' ', '_').replace('/', '-')}_SNPs_pval_0.001.txt"
    if os.path.isfile(f'{gwas_snp_dir}/{gwas_snp_fn}'):
        lip_name = gwas_snp_fn.split('_')[0] # Modified lipid name
        # Get SNPs and dosage
        print(f'\n# Load GWAS SNPs for current lipid: {lip_name}')
        load_dosage_start_time, df_gwas_snp,dosage_all = load_all_dosage(gwas_snp_dir = gwas_snp_dir,
                                                                         gwas_snp_fn = gwas_snp_fn,
                                                                         dosage_dir = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train',
                                                                         dosage_fn = 'species_chr*.vcf.gz.dosage')
        print(f'# - Number of SNPs loaded: {len(df_gwas_snp)}')
        
        print('# Run Elastic net regression')
        # lipid trait, already residuals and looks normal, so no need to INV
        y = df_lipid[lip]
        # y = inverse_normal_transformation(df_lipid[lip])
        # print(y.shape)

        start_time = time.time()
        # Notes from sklearn docs:
        # - l1_ratio is the alpha in R glmnet
        # - alpha is the lambda in R gmlnet
        # Since PrediXcan used glmnet with alpha=0.5,and lambda selected by 10 fold cv,
        # The corresponding parameter in sklearn.ElasticNetCV() are:
        # - l1_ratio=0.5
        # - n_alphas=100, no user supllied selections for alpha, start with n_alphas=10 to save time (#TODO test how long it takes to run a full CV with 100 alphas)
        # - In R glmnet, when nobs > nvars, the default lambda.min.ratio is 0.0001
        # - 10 fold cv
        regr = ElasticNetCV(cv=10,
                            n_alphas=args.n_alphas,
                            random_state=0,
                            n_jobs=32,
                            l1_ratio=0.5) # Default l1 ratio=0.5
        X = dosage_all.T
        regr.fit(X, y)

        end_time = time.time()
        print(f'# - Model fitting finised in {(end_time - start_time):.4f}s')
        
        # Also output predicted values and best R2
        output_fh.write(f"{lip}\t{regr.alpha_}\t{regr.l1_ratio_}\t{regr.score(X, y)}\t{','.join([str(x) for x in regr.coef_])}\n")
        output_fh_lip_pred.write(lip+'\t'+'\t'.join([str(val) for val in regr.predict(X)])+'\n')
        print(f'# Total running time of current lipid: {(end_time - load_dosage_start_time)/60:.4f}m')
        # break
    else:
        print(f'# - Warning: {lip} not found')
    count += 1
    print(f'# #################### {count} lipid processed ####################')
output_fh.close()
output_fh_lip_pred.close()
