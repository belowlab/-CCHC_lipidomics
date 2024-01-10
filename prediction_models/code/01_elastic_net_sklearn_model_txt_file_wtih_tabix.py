# Use tabix to accesss SNPs quickely

# #################################################################
# TODO: need below features
# 1. (DONE) Implement SNP lookups using tabix
# 2. Implement regression types :Elastic net, ridge, lasso and OLS
# #################################################################

from sklearn.linear_model import LassoCV, RidgeCV, ElasticNetCV
import pandas as pd
import numpy as np
import os
import time
import argparse
import subprocess
import gzip
import warnings
import datetime
warnings.filterwarnings(action='ignore')

'''
Example call:
# python 01_elastic_net_sklearn_model_txt_file.py --output lipid_species_l1_0.5_500-599.txt --output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params --lipid_range 500 --range_window 100 --n_alphas 10

python 01_elastic_net_sklearn_model_txt_file.py \
--output output.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/lasso/training/model_params \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/bgziped_dosage \
--dosage_fn species_chr*.vcf.dosage \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3 \
--gwas_snp_fn AC-10:0-_SNPs_pval_0.001.txt \
--lip_name "AC(10:0)" \
--reg_type lasso \
--n_alphas 100 \
--multiallelic False \
--train True

# (Per AlexP) To run in base environment and avoid multi-threading conflict
# Run script as: OMP_NUM_THREADS=1 python my_script.py, for example:

OMP_NUM_THREADS=1 python 01_elastic_net_sklearn_model_txt_file.py \
--output output.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/lasso/training/model_params \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/bgziped_dosage \
--dosage_fn species_chr*.vcf.dosage \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3 \
--gwas_snp_fn AC-10:0-_SNPs_pval_0.001.txt \
--lip_name "AC(10:0)" \
--reg_type lasso \
--n_alphas 100 \
--multiallelic False \
--train True

'''
# ################################## Help functions ##################################

# Print progress par in console
# - progress: current progress (number of SNPs processed)
# - total: total number of SNPs needs to be processed
def progress_bar(progress, total):
    percent = 100 * (progress/total)
    bar = '=' * int(percent) + '-' * int(100 - percent)
    print(f'|{bar}| {percent:.2f}%', end='\r')

def get_doasge(dosage_fn, lst_snp_pos):
    '''
    Use tabix to find dosage of SNPs from a single chromosome
    Param:
     - dosage_fn: name of dosage file to be checked against (single chromosome)
     - lst_snp_pos: start and end position to search for a list of SNP. Such as 1:10000-10001 or chr1:10000-10001.
                Format needs to match chr/pos format in the dosage file
    Return:
     - sample_ids: IDs of genoytpe samples
     - snp_lst: a list of SNPs used in model (ie. SNPs that are found in dosage files)
     - dosage_matrix: dosage of given SNPs as a numpy array. Fill with NA if a SNP is not found
    '''

    # Get genotype sample IDs from dosage file
    with gzip.open(dosage_fn, 'rt') as fh:
        line = fh.readline().strip() # Take sample IDs from header line
        tmp = line.split()
        indx_dosage = tmp.index('FORMAT') + 1 # Get index of sample IDs and dosage values
        sample_ids = tmp[indx_dosage:] # Genotype IDs

    snp_lst = []
    dosage_matrix = []
    count = 0  # Track number of SNPs checked
    for snp_pos in lst_snp_pos:
        tabix_cmd = f'tabix {dosage_fn} {snp_pos}'
        return_vals = subprocess.run(tabix_cmd.split(), capture_output=True, text=True).stdout.strip().split('\n')
        # Drop multiallelic SNPs if args.multiallelic is False
        # If more than one SNPs were found in VCF, ignore this multiallelic site
        if not args.multiallelic:
            if len(return_vals)>1: continue

        try:
            chr_num, pos, snp_id, ref, alt, _, _, _, _, dosage = return_vals[0].split(maxsplit=9)
        except:
            # If multiallelic sites have been already removed from dosage files, tabix call will return empty
            continue

        dosage_matrix.append([float(x) for x in dosage.split()])
        snp_lst.append(snp_id)
        count += 1

        if count%20==0:
            progress_bar(progress=count, total=len(lst_snp_pos))
    progress_bar(progress=1, total=1)
    print(f'\n\t# {count} SNPs processed')

    return sample_ids, snp_lst, np.array(dosage_matrix).reshape(-1, len(sample_ids))

# Load dosage of all filtered SNPs by given p value threshold from GWAS
def load_all_dosage(gwas_snp_fn: str,
                    gwas_snp_dir: str,
                    dosage_dir: str,
                    dosage_fn: str,
                    multiallelic):
    '''
    Get dosage of SNPs from single-chrosmosome dosage files of a given lipid
    Params:
        - gwas_snp_dir: directory to GWAS SNPs (already filtered by p value threshold)
        - gwas_snp_fn: file name of GWAS SNPs (already filtered by p value threshold)
        - dosage_dir: Directory to subsetted dosage file
        - dosage_fn: file name of subset dosage files (by chromosome).
                    Replace chromosome number with '*', such as 'species_chr*.vcf.gz.dosage'
        -multiallelic: drop multiallelic sites if False
    Return:
        - start_time: start time of loading dosage
        - snp_lst: a list of snps loaded
        - dosage_all: A numpy array of dosage. Each row is a SNP, each column is a subject
    '''
    print('# Processing lipid:', args.lip_name)

    # print(f'# Load GWAS SNPs for current lipid')
    df_gwas_snp = pd.read_csv(os.path.join(gwas_snp_dir, gwas_snp_fn), sep='\t').sort_values(by=['CHR', 'POS'])
    # Create regions to lookup usning tabix
    df_gwas_snp['REGION'] = 'chr' + df_gwas_snp['CHR'].astype('str') + ':' + df_gwas_snp['POS'].astype('str') + '-' + df_gwas_snp['POS'].astype('str')

    print('\n# Get dosage of GWAS SNPs to include in regression model')
    print('# - Checking by chromosome:')

    dosage_all = '' # A numpy array to store dosage from all chromosome
    start_time = time.time() # Time execution time
    all_snp_lst = [] # Track what SNPs have been loaded
    for chr_num, df in df_gwas_snp.groupby(by='CHR'):
        # dosage_fn = f'species_chr{chr_num}.vcf.gz.dosage'
        print(f'#  chr{chr_num}')
        sample_ids, snp_lst, dosage_matrix = get_doasge(os.path.join(dosage_dir, dosage_fn.replace('*', str(chr_num))),
                                                        list(df['REGION']))
        # lst_df_dosage.append(pd.DataFrame(data=dosage_matrix, columns=sample_ids, index=df['POS']))
        all_snp_lst += snp_lst
        if len(dosage_all) == 0: # if dosage array is empty
            dosage_all = dosage_matrix
        else:
            dosage_all = np.append(dosage_all, dosage_matrix, axis=0)
        # break
    end_time = time.time()
    print(f'# - Checking finished in {(end_time-start_time):.4f}s')
    print('-' * 50)
    return start_time, all_snp_lst, dosage_all.astype('float64')

# ################################## End of help functions ##################################


# ################################## Process args ##################################
parser = argparse.ArgumentParser(description='Fit regression model of choice (elastic net, ridge, lasso or OLS)')
parser.add_argument('-o', '--output_prefix', type=str,
                           help='Output file to save alpha, l1_ratio and coefficients of chosen model')
parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='.')
parser.add_argument('--dosage_dir', type=str, help='Derictory to dosage files',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train')
parser.add_argument('--dosage_fn', type=str, help='File name format of dosage files. Use * to replace chromosome number',
                    default='species_chr*.vcf.gz.dosage')
parser.add_argument('--gwas_snp_dir', type=str, help='Directory to filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3')
parser.add_argument('--gwas_snp_fn', type=str, help='File name of the filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                    default='AC-10:0-_SNPs_pval_0.001.txt')
parser.add_argument('--lip_name', type=str,
                    help='Name of the lipid to be processed')
parser.add_argument('--n_alphas', type=int, default=100,
                    help='Define how many alphas to test in CV. Dafault is 10. JTI used 100 as defined in R glmnet()')
parser.add_argument('--multiallelic', type=str, default='False',
                    help='If false, multiallelic SNPs will be removed from model fitting')
parser.add_argument('--train', type=str, default='True',
                    help='If true, will not fill with NaN if a SNP is not found. Missing values will cause errors')
parser.add_argument('--reg_type', type=str, default='elastic_net', choices=['elastic_net', 'ridge', 'lasso', 'ols'],
                    help="Type of regression. Choose from: 'elastic_net', 'ridge', 'lasso' and 'ols'")
parser.add_argument('--lipidomis_fn', type=str,
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted',
                    help='Path and name of the lipidomics data to be used in training. Assume values are already transformed or normalized')

# ################################## Arguments sanity checks ##################################
args = parser.parse_args()
args.output_prefix = f"{args.output_prefix}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"
# Check if files exist
for i in range(1, 22):
    if not os.path.isfile(os.path.join(args.dosage_dir, args.dosage_fn.replace('*',str(i)))):
        print('# ERROR: Dosage file not found:', os.path.join(args.dosage_dir, args.dosage_fn))
        # exit()
if not os.path.isfile(os.path.join(args.gwas_snp_dir, args.gwas_snp_fn)):
    print('# ERROR: filtered GWAS result not found:', os.path.join(args.gwas_snp_dir, args.gwas_snp_fn))
    exit()

if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)

if args.multiallelic.upper()[0]=='F' or args.multiallelic=='0':
    args.multiallelic = False
else: # Do not drop multiallelic sites if True
    args.multiallelic = True

if args.train.upper()[0]=='F' or args.train=='0':
    args.train = False
else: # Do not fill missing values with NA
    args.train = True

print('# Run starts:', datetime.datetime.now().strftime('%Y-%m-%d'))
print('# Arguments used:')
for arg in vars(args):
    print(f'# - {arg}:', getattr(args, arg))

# ################################### Load lipidomics data ##################################
print('# Load lipidomic data (lipid species)')
df_lipid = pd.read_csv(args.lipidomis_fn, sep='\t')
print(f"# - data loaded from {args.lipidomis_fn}: shape {df_lipid.shape}")

# Re-order lipidomics data so that sample IDs match the order in genotype file
fn_id_mapping = '/data100t1/home/wanying/CCHC/doc/samples_IDs/202211_merged_RNA_lipid_protein_genotype_mapping_and_availability.txt'
df_id_mapping = pd.read_csv(fn_id_mapping,
                            sep='\t').dropna(subset=['genotype_ID',
                                                     'lipidomic']).drop_duplicates(subset='lipidomic')[['LABID', 'genotype_ID']]

print(f'\n# Load genotype IDs for matching (only need to read the first line of dosage file)')
fn_genotype = os.path.join(args.dosage_dir, args.dosage_fn.replace('*', '22'))
with gzip.open(fn_genotype, 'rt') as fh:
    df_genotype_id = pd.DataFrame(fh.readline().strip().split()[9:], columns=['genotype_ID'])

print(f'# - Organize sample IDs so that their orders match in lipidomics data and dosage file')
df_lipid = df_genotype_id.merge(df_id_mapping.merge(df_lipid.drop_duplicates(subset='Sample ID'),
                                                    left_on='LABID',
                                                    right_on='Sample ID'), on='genotype_ID')
print(f'# - Final processed lipidomics data: {len(df_lipid)}')


# ################################## Load GWAS snps of a given lipid and run regression ##################################
# dosage_all: each row contains dosages of a single SNP across all individuals
# !! Lip species PI(15-MHDA_20:4)\PI(17:0_20:4) is missing
# Save coefficients, alpha and l1 ratios of selected model for each lipid
output_fh = open(os.path.join(args.output_dir, args.output_prefix+f'.{args.reg_type}'), 'w')
output_fh.write('lipid\talpha\tl1_ratio\tbest_r2\tcoefficients\tSNPs\n') # write header line
output_fh_lip_pred = open(os.path.join(args.output_dir, f'{args.output_prefix}.{args.reg_type}.pred'), 'w') # Save predicted values of each lipid using best fitted model
output_fh_lip_pred.write('Lipid'+'\t'+'\t'.join([val for val in df_lipid['Sample ID']])+'\n') # write header line


# #################################################
# Get dosage of a SNP
# Get SNPs and dosage
print(f'\n# Load filtered GWAS SNPs of current lipid: {args.lip_name}')
load_dosage_start_time, snp_lst, dosage_all = load_all_dosage(gwas_snp_dir=args.gwas_snp_dir,
                                                              gwas_snp_fn=args.gwas_snp_fn,
                                                              dosage_dir=args.dosage_dir,
                                                              dosage_fn=args.dosage_fn,
                                                              multiallelic=args.multiallelic)
print(f'# - Number of SNPs loaded: {len(snp_lst)}')

print(f'\n# Run {args.reg_type} regression')
# lipid trait, already residuals and looks normal, so no need to INV
y = df_lipid[args.lip_name]
# y = inverse_normal_transformation(df_lipid[lip])
# print(y.shape)

start_time = time.time()

if args.reg_type == 'elastic_net':
    # Notes from sklearn docs:
    # - l1_ratio is the alpha in R glmnet
    # - alpha is the lambda in R gmlnet
    # Since PrediXcan used glmnet with alpha=0.5,and lambda selected by 10 fold cv,
    # The corresponding parameter in sklearn.ElasticNetCV() are:
    # - l1_ratio=0.5
    # - n_alphas=100, no user supllied selections for alpha, start with n_alphas=10 to save time
    # - In R glmnet, when nobs > nvars, the default lambda.min.ratio is 0.0001
    # - 10 fold cv
    regr = ElasticNetCV(cv=10,
                        n_alphas=args.n_alphas,
                        random_state=0,
                        n_jobs=8,
                        l1_ratio=0.5)  # Default l1 ratio=0.5
elif args.reg_type == 'ridge':
    # regr = ElasticNetCV(cv=10,
    #                     n_alphas=args.n_alphas,
    #                     random_state=0,
    #                     n_jobs=8,
    #                     l1_ratio=0)  # l1 ratio=0 becomes ridge regression
    # TODO: might need this later
    pass
elif args.reg_type == 'lasso':
    regr = ElasticNetCV(cv=10,
                        n_alphas=args.n_alphas,
                        random_state=0,
                        n_jobs=8,
                        l1_ratio=1)  # l1 ratio=1 becomes lasso regression
elif args.reg_type == 'ols':
    # TODO: might need this later
    pass

X = dosage_all.T
regr.fit(X, y)

end_time = time.time()
print(f'# - Model fitting finished in {(end_time - start_time):.4f}s')

# Also output predicted values and best R2
output_fh.write(
    f"{args.lip_name}\t{regr.alpha_}\t{regr.l1_ratio_}\t{regr.score(X, y)}\t{','.join([str(x) for x in regr.coef_])}\t{','.join(snp_lst)}\n")
output_fh_lip_pred.write(args.lip_name + '\t' + '\t'.join([str(val) for val in regr.predict(X)]) + '\n')
print(f'# Total running time of current lipid: {(end_time - load_dosage_start_time) / 60:.4f}m')

print(f'# #################### DONE ####################')

output_fh.close()
output_fh_lip_pred.close()
