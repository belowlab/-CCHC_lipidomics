# This code is to train boosting models on CCHC lipidomics data.
# TODO
# Implement boosting models
# 1. AdaBoost
# 2. Gradient boost
# 3. XGBoost
# Options of base estimator: Linear regression Multilayer perceptron (MLP), Decision tree regressor
# Refer to code ML_03_model_train.py and ML_07_boosting_models_test_run.ipynb
'''
Example call
lip_type=species
lipid="PC(44:5)"
output_prefix=PC-44:5-
python ML_08_boosting_model_train.py \
--output_prefix test_run \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/AdaBoost/${lip_type} \
--dosage_dir_train /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_${lip_type} \
--dosage_fn_train lipid_${lip_type}_chr*.pval_0.001_maf_0.05.vcf.dosage.gz \
--dosage_dir_test /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test/lipid_${lip_type} \
--dosage_fn_test lipid_${lip_type}_chr*.pval_0.001_maf_0.05.test.vcf.dosage.gz \
--gwas_snp_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_${lip_type}_filter_by_pval_1e-07 \
--gwas_snp_fn PC-44:5-_SNPs.pval_1e-07.txt \
--lipid_name ${lipid} \
--trait_fn_train /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.train.txt \
--trait_fn_test /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.test.txt \
--multiallelic False \
--train True \
--n_estimator 10 \
--boost_type Ada
'''

import argparse
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
import logging
import subprocess
import pandas as pd
import numpy as np
from scipy import stats
import os
import datetime
# TODO
# Save model with pickle
# from pickle import dump

# #################### Helper functions ####################
def config_logging(log_fn):
    '''
    Configurate logging. Save messages to an output log file
    :param log_fn:
    :return:
    '''
    logging.getLogger().setLevel(logging.INFO)
    logging.root.handlers = [] # Remove potential handler set up by others
    # logging.info to console and save to log file. Change mode='a' for appending mode
    logging.basicConfig(level=logging.INFO,
                        handlers=[logging.FileHandler(filename=log_fn, mode='a'), logging.StreamHandler()],
                        format='# %(name)s - %(levelname)s - %(message)s')
    # Start logging
    logging.info(__file__)
    logging.info('# ' + '#' * 20 + ' Run started on' + datetime.datetime.now().strftime('%Y-%m-%d') + ' ' +'#' * 20)
    logging.info('')

def parse_arguments():
    '''
    Parse arguments
    :return: Parsed arguments
    '''
    parser = argparse.ArgumentParser(description='Fit regression model of choice (elastic net, ridge, lasso or OLS)')
    parser.add_argument('-o', '--output_prefix', type=str,
                        help='Output file to save alpha, l1_ratio and coefficients of chosen model')
    parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='.')
    parser.add_argument('--dosage_dir_train', type=str, help='Derictory to dosage files of training set',
                        default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_species')
    parser.add_argument('--dosage_fn_train', type=str, help='File name format of dosage files of training set. Use * to replace chromosome number',
                        default='species_chr*.vcf.gz.dosage')
    parser.add_argument('--trait_fn_train', type=str,
                        help='File name of lipidomics data (or other trait) on training set. File must in sample x trait format. Must have a column matches genotype IDs',
                        default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_species_ID_matched.no_dup.residual.train.txt')
    
    # If --trait_fn_test, --dosage_dir_test and --dosage_dir_fn are provided, will apply model on test set and record performance
    parser.add_argument('--dosage_dir_test', type=str, default=None, help='Derictory to dosage files of test set')
    parser.add_argument('--dosage_fn_test', type=str, default=None,
                        help='File name format of dosage files of test set. Use * to replace chromosome number')
    parser.add_argument('--trait_fn_test', type=str,
                        help='File name of lipidomics data (or other trait) of test set. File must in sample x trait format. Must have a column matches genotype IDs',
                        default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_species_ID_matched.no_dup.residual.test.txt')

    parser.add_argument('--gwas_snp_dir', type=str, help='Directory to filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                        default='/data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3')
    parser.add_argument('--gwas_snp_fn', type=str, help='File name of the filtered GWAS SNPs (eg. GWAs SNPs with pval<1e-3)',
                        default='AC-10:0-_SNPs_pval_0.001.txt')
    parser.add_argument('--lipid_name', type=str,
                        help='Name of the lipid to be processed')
    parser.add_argument('--n_estimator', type=int, default=100,
                        help='Define how many estimator to use in the model')
    parser.add_argument('--multiallelic', type=str, default='False',
                        help='If false, multiallelic SNPs will be removed from model fitting')
    parser.add_argument('--train', type=str, default='True',
                        help='If true, will not fill with NaN if a SNP is not found. Missing values will cause errors')
    parser.add_argument('--boost_type', type=str, default='elastic_net', choices=['Ada', 'Gradient', 'XG'],
                        help='Type of boosting. Choose from AdaBoost, gradient and XG')
    return parser.parse_args()

def sanity_checks():
    '''
    Check arguments
    :return:
    '''
    args.output_prefix = f"{args.output_prefix}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"

    # Make all directory to be absolute path
    args.output_dir = os.path.expanduser(args.output_dir)
    args.dosage_dir_train = os.path.expanduser(args.dosage_dir_train)
    args.gwas_snp_dir = os.path.expanduser(args.gwas_snp_dir)
    args.trait_fn_train = os.path.expanduser(args.trait_fn_train)

    msg = '# Check output directory:' # For logging purpose
    if not os.path.isdir(args.output_dir):
        msg += ' Output directory does not exist. Create one at: ' + args.output_dir
        os.makedirs(args.output_dir) # Use makedirs to create folds recursively
    else:
        msg += ' PASS'

    # Start logging
    log_fn = os.path.join(args.output_dir, args.output_prefix+'.log')
    config_logging(log_fn)
    logging.info('# ' + '*' * 20 + ' Sanity checks ' + '*' * 20)
    logging.info(msg)

    # Check if files exist
    logging.info('# Check dosage files of training data')
    all_dosage_files_exist = True
    for i in range(1, 23):
        dosage_fn_train = os.path.join(args.dosage_dir_train, args.dosage_fn_train.replace('*', str(i)))
        if not os.path.isfile(dosage_fn_train):
            logging.error('# - CHR%s dosage file not found: %s' % (i, dosage_fn_train))
            all_dosage_files_exist = False
        else:
            logging.info('# - CHR%s: PASS' % i)
    if not all_dosage_files_exist:
        logging.error('# - Missing dosage files. Exit')
        exit()

    if (args.dosage_dir_test is not None) and (args.dosage_fn_test is not None):
        logging.info('# Check dosage files of test data')
        all_dosage_files_exist = True
        dosage_fn_test = os.path.join(args.dosage_dir_test, args.dosage_fn_test.replace('*', str(i)))
        if not os.path.isfile(dosage_fn_test):
            logging.error('# - CHR%s dosage file not found: %s' % (i, dosage_fn_test))
            all_dosage_files_exist = False
        else:
            logging.info('# - CHR%s: PASS' % i)
        if not all_dosage_files_exist:
            logging.error('# - Missing dosage files. Exit')
            exit()
        
    logging.info('# Check (filtered) GWAS SNP file: ')
    if not os.path.isfile(os.path.join(args.gwas_snp_dir, args.gwas_snp_fn)):
        logging.error('# - filtered GWAS result not found: ' + os.path.join(args.gwas_snp_dir, args.gwas_snp_fn))
        exit()
    else:
        logging.info('# - PASS')

    logging.info('# Check trait file (residuals of lipidomic measures):')
    if not os.path.isfile(args.trait_fn_train):
        logging.error('# - trait file not found: ' + args.trait_fn_train)
        exit()
    else:
        logging.info('# - PASS')

    if args.multiallelic.upper()[0] == 'F' or args.multiallelic == '0':
        args.multiallelic = False
    else:  # Do not drop multiallelic sites if True
        args.multiallelic = True

    if args.train.upper()[0] == 'F' or args.train == '0':
        args.train = False
    else:  # Do not fill missing values with NA # Seems unnecessary?
        args.train = True
    logging.info('# Arguments used:')
    for arg in vars(args):
        logging.info('# - %s: %s' % (arg, getattr(args, arg)))
    return msg

def load_dosage(snp_dir, snp_fn, dosage_dir, dosage_fn):
    '''
    Load dosage data of a given lipid trait
    :param snp_dir:
    :param snp_fn:
    :param dosage_dir:
    :param dosage_fn:
    :return:
    '''
    logging.info('# Load filtered SNPs')
    df_snps = pd.read_csv(os.path.join(snp_dir, snp_fn), sep='\t')

    logging.info('# Get SNP dosage from CHR1-22')
    snps_all, dosage_all, lst_sample_ids = [], [], [] # lists of SNPs and dosage found and used in model, sample IDs
    total_num_snps = len(df_snps)
    c_all_loaded = 0 # track number of all loaded snps
    for chr_num, df in df_snps.groupby(by='CHR'):
        logging.info('# - Get dosage from CHR%s'%chr_num)
        # Create a list of positions to search by tabix
        lst_pos = 'chr' + df['CHR'].astype('str') + ':' + df['POS'].astype('str') + '-' + df['POS'].astype('str')
        total, c = len(df), 0 # Number of SNPs to be found and actually loaded  on current chromosome
        count_missing = 0  # Count number of missing or skipped sites
        for pos in lst_pos:
            cmd = f'tabix {os.path.join(dosage_dir, dosage_fn.replace("*", str(chr_num)))} {pos} | cut -f 1-5,10-'
            dosage = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout.strip()
            if dosage == '':
                count_missing += 1
                continue # If the position is not found, tabix call returns empty string
            if not args.multiallelic: # Ignore multiallelic sites when args.multiallelic is False
                if len(dosage.split('\n')) > 1:
                    count_missing += 1
                    continue
            tmp_lst = dosage.split('\n')[0].split()
            snps_all.append(tmp_lst[2])
            dosage_all.append(tmp_lst[5:])
            c += 1
            if c == 1: # Get sample IDs from header by using -H flag with tabix
                cmd = f'tabix -H {os.path.join(dosage_dir, dosage_fn.replace("*", str(chr_num)))} {pos} | cut -f 10-'
                sample_ids = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout
                lst_sample_ids = sample_ids.split('\n')[-2].split('\t')
        logging.info('#    - CHR%s: SNPs loaded %s/%s; Multiallelic or missing SNPs: %s' % (chr_num, c, total, count_missing))
        c_all_loaded += c
    df_dosage_all = pd.DataFrame(data=np.array(dosage_all), columns=lst_sample_ids, index=snps_all)
    logging.info('# - Total number of SNPs loaded: %s' % c)
    return df_dosage_all.T.reset_index().rename(columns={'index':'genotype_ID'})


# TODO
# fix below code
def load_trait_values(fn_trait):
    '''
    Load trait values (lipidomic residuals)
    :param fn_trait: file name
    :return: a DataFrame of trait values, a dictionary matches lipid and lipid names
    '''
    # create a dictionary for modified lipid name matching
    # Lipid and lipid name are not the same. I did a few replacement to avoid special characters in file name
    # TODO: lipid trait list is hard coded, might need to change, so it is more flexible
    # fn_list_trait = f'/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_{lipid_type}.list'
    # df_trait_name_matching = pd.read_csv(fn_list_trait, sep='\t', header=None).rename(columns={0:'Lipid'})
    # df_trait_name_matching['Lipid_name'] = df_trait_name_matching['Lipid'].apply(lambda x: x.replace('\\', '-').replace('/', '-').replace('(','-').replace(')','-').replace(' ', '_'))
    # dict_trait_name_matching = df_trait_name_matching.set_index(keys='Lipid_name').to_dict()['Lipid']
    
    logging.info('# Load trait values of all lipids')
    # trait_dir = '/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait'
    # trait_fn = f'lipid_{lipid_type}_ID_matched.no_dup.residual.{data_type}.txt'
    df_trait = pd.read_csv(fn_trait, sep='\t')
    return df_trait

# #################### End of helper functions ####################

args = parse_arguments()
sanity_checks()

logging.info('# ' + '*' * 20 + ' Get dosage ' + '*' * 20)
# Load dosage
df_dosage = load_dosage(snp_dir=args.gwas_snp_dir,
                        snp_fn=args.gwas_snp_fn,
                        dosage_dir=args.dosage_dir_train,
                        dosage_fn=args.dosage_fn_train)

# #################### Load trait ####################
logging.info('# ' + '*' * 20 + ' Get trait values ' + '*' * 20)
df_trait = load_trait_values(fn_trait=args.trait_fn_train)

# Reorder samples so that dosage and trait dataframe match
logging.info('# - Order samples in dosage and trait so that they match')
assert 'genotype_ID' in df_dosage.columns
assert 'genotype_ID' in df_trait.columns
if len(df_trait)>len(df_dosage):
    # Take sample IDs in the smaller dataframe to avoid NA
    samples_index=df_dosage['genotype_ID']
    df_trait = df_trait.set_index(keys='genotype_ID').reindex(samples_index).reset_index()
else:
    samples_index=df_trait['genotype_ID']
    df_dosage = df_dosage.set_index(keys='genotype_ID').reindex(samples_index).reset_index()

# #################### Train model ####################
logging.info('# ' + '*' * 20 + ' Model training ' + '*' * 20)
X = df_dosage.iloc[:, 1:].values
y = df_trait[args.lipid_name]
assert X.shape[0] == y.shape[0]
logging.info('# %s samples were used in model training' % y.shape[0])

n_estimators = args.n_estimator
regr = AdaBoostRegressor(random_state=0, n_estimators=n_estimators)
regr.fit(X, y)
logging.info('# N estimators=%s' % args.n_estimator)
logging.info('# Model fitting r2: %.4f' % regr.score(X, y))

# #################### Apply model on test set ####################
logging.info('# ' + '*' * 20 + ' Apply model on test set ' + '*' * 20)
logging.info('# Get dosage')
# Load dosage
df_dosage_test = load_dosage(snp_dir=args.gwas_snp_dir,
                             snp_fn=args.gwas_snp_fn,
                             dosage_dir=args.dosage_dir_test,
                             dosage_fn=args.dosage_fn_test)

# Load trait of test set
logging.info('# - Get trait values')
df_trait = load_trait_values(fn_trait=args.trait_fn_test)

# Reorder samples so that dosage and trait dataframe match
logging.info('# - Order samples in dosage and trait so that they match')
assert 'genotype_ID' in df_dosage_test.columns
assert 'genotype_ID' in df_trait.columns
if len(df_trait)>len(df_dosage_test):
    # Take sample IDs in the smaller dataframe to avoid NA
    samples_index=df_dosage_test['genotype_ID']
    df_trait = df_trait.set_index(keys='genotype_ID').reindex(samples_index).reset_index()
else:
    samples_index=df_trait['genotype_ID']
    df_dosage_test = df_dosage_test.set_index(keys='genotype_ID').reindex(samples_index).reset_index()
    
X_test = df_dosage_test.iloc[:, 1:].values
y_test = df_trait[args.lipid_name]

print(df_dosage_test)
print(df_dosage_test.shape)
print(y_test.shape)


logging.info('# %s samples were used in testing' % y_test.shape[0])
assert X_test.shape[0] == y_test.shape[0]
logging.info('# Pearson r2: %.4f' % stats.pearsonr(y_test, regr.predict(X_test))[0]**2)
