# Modified from ~/CCHC/lipidomics/prediction_models/code/01_elastic_net_sklearn_model_txt_file_with_tabix.py
# Extract SNP dosage and train model of a given lipid trait

# #################################################################
# TODO: need below features
# 1. (DONE) Implement SNP lookups using tabix
# 2. Implement regression types: Elastic net, ridge, lasso and OLS
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
from scipy import stats

'''
Example call:
# (Per AlexP) To run in base environment and avoid multi-threading conflict
# Run script as: OMP_NUM_THREADS=1 python my_script.py, for example:

python ML_03_model_train.py \
--output_prefix test_run \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/class \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_class \
--dosage_fn lipid_class_chr*.pval_0.001_maf_0.05.vcf.dosage.gz \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_class_filter_by_pval_1e-05 \
--gwas_snp_fn AC-OH_SNPs.pval_1e-05.txt \
--lip_name "AC(10:0)" \
--n_alphas 100 \
--multiallelic False \
--train True \
--reg_type elastic_net

lip_type=species
lipid="PC(44:5)"
python ML_03_model_train.py \
--output_prefix test_run \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/${lip_type} \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_${lip_type} \
--dosage_fn lipid_${lip_type}_chr*.pval_0.001_maf_0.05.vcf.dosage.gz \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_${lip_type}_filter_by_pval_1e-05 \
--gwas_snp_fn PC-44:5-_SNPs.pval_1e-05.txt \
--lip_name ${lipid} \
--trait_fn /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.train.txt \
--n_alphas 100 \
--multiallelic False \
--train True \
--reg_type elastic_net


OMP_NUM_THREADS=1 python 01_elastic_net_sklearn_model_txt_file.py \
--output output.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/lasso/training/model_params \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/subset_vcfs/train/bgziped_dosage \
--dosage_fn species_chr*.vcf.dosage \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3 \
--gwas_snp_fn AC-10:0-_SNPs_pval_0.001.txt \
--lip_name "AC(10:0)" \
--trait_fn /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_species_ID_matched.no_dup.residual.train.txt \
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

def load_trait(trait_fn, genotype_ids):
    '''
    Load dependent variable Y (ie. residuals of lipidomics)
    :param trait_fn: file name (of lipidomics data or other trait).
                     File must in sample x trait format.
                     Must have a column matches genotype IDs
           genotype_ids: array-like genotype IDs for re-ordering
    :return: a dataframe of traits to be used for model training
    '''
    if trait_fn.endswith('.csv'):
        df = pd.read_csv(trait_fn)
    else:
        df = pd.read_csv(trait_fn, sep='\t')
        
    # Re-order lipidomics data so that sample IDs match the order in genotype file
    # df.set_index(keys='genotyoe_IDs', inplace=True)
    return df.set_index(keys='genotype_ID').reindex(index=genotype_ids).reset_index()

def get_dosage(dosage_fn, lst_snp_pos, chr_num):
    '''
    Use tabix to find dosage of SNPs from a single chromosome.
    Need to bgzip and tabix the dosage files first
    Param:
     - dosage_fn: name of dosage file to be checked against (single chromosome).
                  Use filtered dosage files to speed up. Tabix index dosage file before model training.
                  Dosage file must have one header line of sample IDs.
     - lst_snp_pos: start and end position (can be the same) to search for a list of SNPs. Such as 1:10000-10000 or chr1:10000-10000.
                    Format needs to match chr/pos format in the dosage file
     - chr_num: chromosome number, for console message
    Return:
     - sample_ids: IDs of genotype samples
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
    total = len(lst_snp_pos)
    count_multiallelic = 0 # Count number of multiallelic sites skipped
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
            count_multiallelic += 1
            continue

        dosage_matrix.append([float(x) for x in dosage.split()])
        snp_lst.append(snp_id)
        count += 1
        print(f'\r# - {chr_num}: {count}/{total} SNPs loaded; {count_multiallelic} not found (multiallelic sites)  ',
              end='', flush=True)
    print(f'\r# - {chr_num}: {count}/{total} SNPs loaded; {count_multiallelic} not found (multiallelic sites)  ', flush=True)

    return sample_ids, snp_lst, np.array(dosage_matrix).reshape(-1, len(sample_ids))

# Load dosage of all filtered SNPs from all chromosomes by given p value (or MAF) threshold from GWAS
def load_all_dosage(gwas_snp_fn: str,
                    gwas_snp_dir: str,
                    dosage_dir: str,
                    dosage_fn: str):
    '''
    Get dosage of SNPs from single-chrosmosome dosage files of a given lipid
    Params:
        - gwas_snp_dir: directory to GWAS SNPs (already filtered by p value, MAF threshold, etc)
        - gwas_snp_fn: file name of GWAS SNPs (already filtered by p value, MAF threshold, etc.)
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
    # Create regions to lookup using tabix
    df_gwas_snp['REGION'] = 'chr' + df_gwas_snp['CHR'].astype('str') + ':' + df_gwas_snp['POS'].astype('str') + '-' + df_gwas_snp['POS'].astype('str')

    print('\n# Get dosage of GWAS SNPs to include in regression model')
    print('# Checking by chromosome:')

    dosage_all = '' # A numpy array to store dosage from all chromosome
    start_time = time.time() # Time execution time
    all_snp_lst = [] # Track what SNPs have been loaded
    for chr_num, df in df_gwas_snp.groupby(by='CHR'):
        # dosage_fn = f'species_chr{chr_num}.vcf.gz.dosage'
        # print(f'# - chr{chr_num}: ', end='')
        sample_ids, snp_lst, dosage_matrix = get_dosage(os.path.join(dosage_dir, dosage_fn.replace('*', str(chr_num))),
                                                        list(df['REGION']), chr_num)
        # lst_df_dosage.append(pd.DataFrame(data=dosage_matrix, columns=sample_ids, index=df['POS']))
        all_snp_lst += snp_lst
        if len(dosage_all) == 0: # if dosage array is empty
            dosage_all = dosage_matrix
        else:
            dosage_all = np.append(dosage_all, dosage_matrix, axis=0)
    return start_time, all_snp_lst, dosage_all.astype('float64')

# ################################## End of help functions ##################################


# ################################## Process args ##################################
parser = argparse.ArgumentParser(description='Fit regression model of choice (elastic net, ridge, lasso or OLS)')
parser.add_argument('-o', '--output_prefix', type=str,
                           help='Output file to save alpha, l1_ratio and coefficients of chosen model')
parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='.')
parser.add_argument('--dosage_dir', type=str, help='Derictory to dosage files',
                    default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_species')
parser.add_argument('--dosage_fn', type=str, help='File name format of dosage files. Use * to replace chromosome number',
                    default='species_chr*.vcf.gz.dosage')
parser.add_argument('--gwas_snp_dir', type=str, help='Directory to filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3')
parser.add_argument('--gwas_snp_fn', type=str, help='File name of the filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                    default='AC-10:0-_SNPs_pval_0.001.txt')
parser.add_argument('--lip_name', type=str,
                    help='Name of the lipid to be processed')
parser.add_argument('--trait_fn', type=str, help='File name of lipidomics data (or other trait). File must in sample x trait format. Must have a column matches genotype IDs',
                    default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_species_ID_matched.no_dup.residual.train.txt')
parser.add_argument('--n_alphas', type=int, default=100,
                    help='Define how many alphas to test in CV. Default is 10. JTI used 100 as defined in R glmnet()')
parser.add_argument('--multiallelic', type=str, default='False',
                    help='If false, multiallelic SNPs will be removed from model fitting')
parser.add_argument('--train', type=str, default='True',
                    help='If true, will not fill with NaN if a SNP is not found. Missing values will cause errors')
parser.add_argument('--reg_type', type=str, default='elastic_net', choices=['elastic_net', 'ridge', 'lasso', 'ols'],
                    help="Type of regression. Choose from: 'elastic_net', 'ridge', 'lasso' and 'ols'")

# ################################## Arguments sanity checks ##################################
print('# Run starts:', datetime.datetime.now().strftime('%Y-%m-%d'))
print('\n')

print('#', '#'*40, 'Sanity checks', '#'*40)
args = parser.parse_args()
args.output_prefix = f"{args.output_prefix}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"

# Make all directory to absolute path
args.output_dir = os.path.expanduser(args.output_dir)
args.dosage_dir = os.path.expanduser(args.dosage_dir)
args.gwas_snp_dir = os.path.expanduser(args.gwas_snp_dir)
args.trait_fn = os.path.expanduser(args.trait_fn)

# Check if files exist
print('# - Check dosage files')
all_dosage_files_exist = True
for i in range(1, 22):
    if not os.path.isfile(os.path.join(args.dosage_dir, args.dosage_fn.replace('*', str(i)))):
        print(f'#\t - chr{i}: Dosage file not found:', os.path.join(args.dosage_dir, args.dosage_fn.replace('*', str(i))))
        all_dosage_files_exist = False
    else: print(f'#\t - chr{i}: PASS')
if not all_dosage_files_exist:
    print('# - Missing dosage files. Exit')
    exit()

print('# - Check (filtered) GWAS SNP file: ', end='')
if not os.path.isfile(os.path.join(args.gwas_snp_dir, args.gwas_snp_fn)):
    print('filtered GWAS result not found:', os.path.join(args.gwas_snp_dir, args.gwas_snp_fn))
    exit()
else:
    print('PASS')

print('# - Check trait file (residuals of lipidomic measures): ', end='')
if not os.path.isfile(args.trait_fn):
    print('trait file not found:', args.trait_fn)
    exit()
else:
    print('PASS')

print('# - Check output directory: ', end='')
if not os.path.isdir(args.output_dir):
    print('\n#\t - Output directory does not exist. Create one at:', args.output_dir)
    os.mkdir(args.output_dir)
else:
    print('PASS')

if args.multiallelic.upper()[0]=='F' or args.multiallelic=='0':
    args.multiallelic = False
else: # Do not drop multiallelic sites if True
    args.multiallelic = True

if args.train.upper()[0]=='F' or args.train=='0':
    args.train = False
else: # Do not fill missing values with NA
    args.train = True

print('\n# Arguments used:')
for arg in vars(args):
    print(f'# - {arg}:', getattr(args, arg))

# ################################### Load lipidomics data ##################################
print('\n#', '#'*40, 'Load lipidomics data', '#'*40)
# Re-order lipidomics data so that sample IDs match the order in genotype file
# print(f'\n# Load genotype IDs for matching (only need to read the first line of dosage file)')
fn_genotype = os.path.join(args.dosage_dir, args.dosage_fn.replace('*', '22'))
with gzip.open(fn_genotype, 'rt') as fh:
    genotype_ids = fh.readline().strip().split()[9:]
df_lipid = load_trait(args.trait_fn, genotype_ids)
print(f"# - data loaded from {args.trait_fn}: shape {df_lipid.shape}")
print(f'# - Final processed lipidomics data: {len(df_lipid)}')

# ################################## Load GWAS snps of a given lipid and run regression ##################################
# dosage_all: each row contains dosages of a single SNP across all individuals
# Save coefficients, alpha and l1 ratios of selected model for each lipid
output_fh = open(os.path.join(args.output_dir, args.output_prefix+f'.{args.reg_type}'), 'w')
output_fh.write('lipid\talpha\tl1_ratio\tbest_r2\tcoefficients\tSNPs\n') # write header line
output_fh_lip_pred = open(os.path.join(args.output_dir, f'{args.output_prefix}.{args.reg_type}.pred'), 'w') # Save predicted values of each lipid using best fitted model
output_fh_lip_pred.write('Lipid'+'\t'+'\t'.join([val for val in df_lipid['genotype_ID']])+'\n') # write header line

# Get SNPs and dosage
print('\n#', '#'*40, f'Load filtered GWAS SNPs of current lipid: {args.lip_name}', '#'*40)
start_time = time.time()
load_dosage_start_time, snp_lst, dosage_all = load_all_dosage(gwas_snp_dir=args.gwas_snp_dir,
                                                              gwas_snp_fn=args.gwas_snp_fn,
                                                              dosage_dir=args.dosage_dir,
                                                              dosage_fn=args.dosage_fn)
end_time = time.time()
print('-' * 50)
print(f'# - Number of SNPs loaded: {len(snp_lst)}')
print(f'# - Variants lookup time: {(end_time - start_time):.4f}s')

start_time = time.time()
print('\n#', '#'*40, 'Model training', '#'*40)
print(f'# Run {args.reg_type} regression')

y = df_lipid[args.lip_name].values
# Drop missing values in genotype file
indx = np.argwhere(np.isnan(y))
if indx != 0: # If there are missing values
    print(f'# - Drop missing values in trait: N={len(indx)}')
y = np.delete(y, indx)

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
    # Remember to save Intercept!!!!
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
if len(indx) != 0: # Drop rows in dosage according to missing values in trait
    X = np.delete(X, indx, axis=0)
    print(f'# - Drop missing values in dosage: N={len(indx)}')
regr.fit(X, y)

end_time = time.time()
print(f'# - Model fitting finished in {(end_time - start_time):.4f}s')

# Also output predicted values and best R2
# Need to save the intercept!!!
output_fh.write(
    f"{args.lip_name}\t{regr.alpha_}\t{regr.l1_ratio_}\t{regr.score(X, y)}\t{','.join([str(x) for x in regr.coef_])}\t{','.join(snp_lst)}\n")
output_fh_lip_pred.write(args.lip_name + '\t' + '\t'.join([str(val) for val in regr.predict(X)]) + '\n')
print(f'# Total running time of current lipid: {(end_time - load_dosage_start_time) / 60:.4f}m')

print('\n#', '#'*40, 'DONE', '#'*40)
# print('#### INTERCEPT', regr.intercept_)
output_fh.close()
output_fh_lip_pred.close()



y_pred = regr.predict(X)
print('# Pearson r=%s, pval=%s' % stats.pearsonr(y, y_pred))

# #############
print('#', '-'*50)
print('# y true:', len(y))
print(','.join([str(val) for val in y]))
print('# y pred:', len(y_pred))
print(','.join([str(val) for val in y_pred]))
print()