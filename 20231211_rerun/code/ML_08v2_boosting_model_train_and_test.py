# Modified from ML_08_boosting_model_train.py
# Train and test boosting models. Use 10-fold CV

'''
Example call
lip_type=species
lipid="AC(16:1)-OH"
output_prefix=AC-16:1--OH
python ML_08v2_boosting_model_train_and_test.py \
--output_prefix ${output_prefix} \
--output_dir ./ \
--dosage_dir_train /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_${lip_type} \
--dosage_fn_train lipid_${lip_type}_chr*.pval_0.001_maf_0.05.vcf.dosage.gz \
--dosage_dir_test /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test/lipid_${lip_type} \
--dosage_fn_test lipid_${lip_type}_chr*.pval_0.001_maf_0.05.test.vcf.dosage.gz \
--gwas_snp_fn /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_species_filter_by_pval_1e-07/AC-16:1--OH_SNPs.pval_1e-07.txt \
--lipid_name ${lipid} \
--trait_fn_train /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.train.txt \
--trait_fn_test /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.test.txt \
--multiallelic False \
--n_estimator 20,30,50,75,100,150,200,250,300 \
--learning_rate 0.002,0.004,0.006,0.008,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.5 \
--model Gradient
'''

from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, ElasticNet, Lasso
from sklearn.metrics import r2_score
from sklearn.experimental import enable_halving_search_cv # Required for HalvingGridSearch
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.model_selection import GridSearchCV
import argparse
import logging
import subprocess
from scipy import stats
import pandas as pd
import numpy as np
import os
import datetime

import sys
sys.path.insert(0, '/data100t1/home/wanying/CCHC/lipidomics/code/utils') # Search util first to use my functions
# Import functions from my module
import load_dosage
import load_trait_values
import set_logging

# Save model with pickle also works, but not the preferred way
# from pickle import dump
from joblib import dump, load # To save trained model

# ###################### Helper function ######################
def parse_arguments():
    '''
    Parse arguments
    :return: Parsed arguments
    '''
    parser = argparse.ArgumentParser(description='Fit ensemble (boosting) model of choice (AdaBoost, GradientBoosting)')
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
                        default='')

    parser.add_argument('--gwas_snp_fn', type=str, help='Path and file name of the filtered GWAS SNPs (eg. SNPs with GWAS pval<1e-3)',
                        default=None)
    parser.add_argument('--lipid_name', type=str,
                        help='Name of the single lipid to be processed. For example: PC(44:5). Depend onthe name used in trait file')
    parser.add_argument('--n_estimator', type=str, default='100',
                        help='Define how many estimator to use in the model. Can be a single value or a list of values separated by comma. Eg. 50, or 20,50,100')
    parser.add_argument('--learning_rate', type=str, default='1',
                        help='Define learning rate of the model the model. Can be a single value (between 0 to 1) or a list of values separated by comma. Eg. 1 or 0.1,0.3,0.5')
    parser.add_argument('--multiallelic', type=str, default='False',
                        help='If false, multiallelic SNPs will be removed from model fitting')
    parser.add_argument('--model', type=str, default='Ada', choices=['Ada', 'Gradient', 'XG'],
                        help='Type of boosting. Choose from AdaBoost, gradient and XG (have not implemented XGboosting yet)')
    return parser.parse_args()

def sanity_checks():
    '''
    Check arguments
    :return: Exit script if anything went wrong
    '''
    args.output_prefix = f"{args.output_prefix}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"

    # Make all directory to be absolute path
    args.output_dir = os.path.expanduser(args.output_dir)
    args.dosage_dir_train = os.path.expanduser(args.dosage_dir_train)
    args.trait_fn_train = os.path.expanduser(args.trait_fn_train)
    args.gwas_snp_fn = os.path.expanduser(args.gwas_snp_fn)

    msg = '# Check output directory:' # For logging purpose
    if not os.path.isdir(args.output_dir):
        msg += ' Output directory does not exist. Create one at: ' + args.output_dir
        os.makedirs(args.output_dir) # Use makedirs to create folds recursively
    else:
        msg += ' PASS'

    # Start logging
    log_fn = os.path.join(args.output_dir, args.output_prefix+'.log')
    set_logging.set_logging(log_fn)
    logging.info('# ' + '*' * 20 + ' Sanity checks ' + '*' * 20)
    logging.info(msg)

    # Check if files exist
    logging.info('# Check dosage files of training data')
    all_dosage_files_exist = True
    missing_chr = '' # Track missing chromosome
    for i in range(1, 23):
        dosage_fn_train = os.path.join(args.dosage_dir_train, args.dosage_fn_train.replace('*', str(i)))
        if not os.path.isfile(dosage_fn_train):
            logging.error('# - CHR%s dosage file not found: %s' % (i, dosage_fn_train))
            all_dosage_files_exist = False
            missing_chr += ' '+str(i)
    if not all_dosage_files_exist:
        logging.error('# - Missing dosage files of %s. Exit' % missing_chr)
        exit()
    logging.info('# - PASS')
    
    logging.info('# Check trait file of training set(residuals of lipidomic measures):')
    if not os.path.isfile(args.trait_fn_train):
        logging.error('# - trait file not found: ' + args.trait_fn_train)
        exit()
    else: logging.info('# - PASS')
    
    # Apply model on test set if dosage files of test set are provided
    if (args.dosage_dir_test is not None) and (args.dosage_fn_test is not None):
        logging.info('# Check dosage files of test data')
        all_dosage_files_exist = True
        missing_chr += ''
        for i in range(1, 23):
            dosage_fn_test = os.path.join(args.dosage_dir_test, args.dosage_fn_test.replace('*', str(i)))
            if not os.path.isfile(dosage_fn_test):
                logging.error('# - CHR%s dosage file not found: %s' % (i, dosage_fn_test))
                all_dosage_files_exist = False
                missing_chr += ' '+str(i)
        if not all_dosage_files_exist:
            logging.error('# - Missing dosage files of test set. Exit')
            exit()
        logging.info('# - PASS')
        logging.info('# Check trait files of test data')
        args.trait_fn_test = os.path.expanduser(args.trait_fn_test)
        if not os.path.isfile(args.trait_fn_test):
            logging.error('# - Trait files of test set not found. Exit')
            exit()
        else: logging.info('# - PASS')
        
    logging.info('# Check (filtered) GWAS SNP file: ')
    if not os.path.isfile(args.gwas_snp_fn):
        logging.error('# - filtered GWAS result not found: ' + args.gwas_snp_fn)
        exit()
    else: logging.info('# - PASS')

    if args.multiallelic.upper()[0] == 'F' or args.multiallelic == '0':
        args.multiallelic = False
    else:  # Do not drop multiallelic sites if True
        args.multiallelic = True
        
    # Process n_estimator and learning_rate
    logging.info('# Process n_estimator and learning_rate')
    args.n_estimator = [int(x) for x in args.n_estimator.split(',')]
    args.learning_rate = [float(x) for x in args.learning_rate.split(',')]
    for lr in args.learning_rate: # learning rate should between 0 and 1
        if lr>1 or lr<0:
            logging.error('# - Invalid learning rate')
            exit()
    logging.info('# - PASS')
        
    logging.info('# Arguments used:')
    for arg in vars(args):
        if args in ['n_estimator', 'learning_rate']:
            msg = ', '.join(getattr(args, arg))
            logging.info('# - %s: %s' % (arg, msg))
        else:
            logging.info('# - %s: %s' % (arg, getattr(args, arg)))
    return msg

def prepare_data(snp_fn, dosage_dir, dosage_fn, trait_fn):
    '''
    Load dosafe and trait. Get X and y for training or prediction
    :param snp_fn: file name of filtered GWAS snps
    :param dosage_dir:
    :param dosage_fn:
    :param trait_fn: file name of lipid trait
    
    :return:
    df_dosage, df_lipid_trait: re-ordered dataframes
    X: Features
    y: outcome
    '''
    # Load dosage
    df_dosage = load_dosage.load_dosage(snp_dir='',
                                        snp_fn=snp_fn,
                                        dosage_dir=dosage_dir,
                                        dosage_fn=dosage_fn)
    # Load lipid trait
    df_lipid_trait = load_trait_values.load_trait_values(trait_dir='',
                                                         trait_fn=trait_fn)
    logging.info('# - Reorder samples, remove sample with missing value')
    # Reorder dosage and trait dataframes
    # Take the shorter one between dosage and lipid measure to avoid missing values
    if len(df_dosage) < len(df_lipid_trait):
        df_lipid_trait = df_lipid_trait.set_index(keys='genotype_ID').reindex(df_dosage['genotype_ID']).reset_index()
    else:
        df_dosage = df_dosage.set_index(keys='genotype_ID').reindex(df_lipid_trait['genotype_ID']).reset_index()

    # Create X, y for model training
    X = df_dosage.iloc[:, 1:].values
    y = df_lipid_trait[args.lipid_name]
    return df_dosage, df_lipid_trait, X, y
# #################### End of helper functions ####################

start = datetime.datetime.now()
args = parse_arguments()
sanity_checks()

logging.info('')
logging.info('# Load data')
df_dosage, df_lipid_trait, X, y = prepare_data(snp_fn=args.gwas_snp_fn,
                                               dosage_dir=args.dosage_dir_train,
                                               dosage_fn=args.dosage_fn_train,
                                               trait_fn=args.trait_fn_train)
logging.info('# - Shape of X, y in training: %s, %s' % (X.shape, y.shape))

logging.info('')
logging.info('# Start model train')
logging.info('# - Use default estimator: DT regressor')
if args.model == 'Ada':
    logging.info('# - Use AdaBoost')
    regr = AdaBoostRegressor(random_state=0)
elif args.model == 'Gradient':
    logging.info('# - Use gradient boosting')
    regr = GradientBoostingRegressor(random_state=0)
else:
    logging.info('# - Method has not been implemented yet. Exit')
    exit()
    
param_grid = {'n_estimators':args.n_estimator,
              'learning_rate':args.learning_rate}

# Use HalvingGridSearchCV as it is faster than GridSearchCV
logging.info('# - Use HalvingGridSearchCV, which is faster than GridSearchCV')
search = HalvingGridSearchCV(estimator=regr, param_grid=param_grid, cv=10, n_jobs=16, random_state=1).fit(X, y)
# Summarize the best score and configuration
logging.info('# Model fitting metrics')
logging.info("# - Best: socre r2=%f using %s" % (search.best_score_, search.best_params_))

# Summarize all scores that were evaluated
means = search.cv_results_['mean_test_score']
stds = search.cv_results_['std_test_score']
params = search.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    logging.info("# - Coefficient of determination score (validation) r2=%f (std=%f) with: %r" % (mean, stdev, param))

# #################### Save model and predicted values ####################
logging.info('# ' + '*' * 20 + ' Save model and predicted values ' + '*' * 20)
cv_results_fn = os.path.join(args.output_dir, args.output_prefix)+'.cv_results'
logging.info('# Model training metrics saved to %s' % cv_results_fn)
# Save training metrics for future reference
df_halving = pd.DataFrame(search.cv_results_) # Output some metrics
df_halving.to_csv(cv_results_fn, sep='\t', index=False)

fn_model = os.path.join(args.output_dir, args.output_prefix+'.joblib')
dump(search, filename=fn_model)
logging.info('# Model saved to: %s' % fn_model)

# Apply on test set if test dosage files are provided
if (args.dosage_dir_test is not None) and (args.dosage_fn_test is not None):
    logging.info('')
    logging.info('#' + '*' * 20 + ' Apply model on test set ' + '*' * 20)

    df_dosage_test, df_lipid_trait_test, X_test, y_test = prepare_data(snp_fn=args.gwas_snp_fn,
                                                                       dosage_dir=args.dosage_dir_test,
                                                                       dosage_fn=args.dosage_fn_test,
                                                                       trait_fn=args.trait_fn_test)
    logging.info('# - Shape of X, y in testing: %s, %s' % (X_test.shape, y_test.shape))
    y_pred_test = search.predict(X_test)
    fn_pred = os.path.join(args.output_dir, args.output_prefix+'.pred')
    logging.info('# Coefficient of determination score (test) r2=%.4f')
    logging.info('Pearson r2 on test: %.4f', stats.pearsonr(y_test, y_pred_test)[0]**2)
    
    y_pred_test = [str(x) for x in y_pred_test]
    with open(fn_pred, 'w') as fh:
        fh.write('\t'.join(df_lipid_trait_test['genotype_ID'])+'\n')
        fh.write('\t'.join(y_pred_test))
    logging.info('# Predicted values saved to: %s' % fn_pred)