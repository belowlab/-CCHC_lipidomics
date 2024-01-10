# This code is to train boosting models on CCHC lipidomics data.
# TODO
# Implement boosting models
# 1. AdaBoost
# 2. Gradient boost
# 3. XGBoost
# Options of base estimator: Linear regression Multilayer perceptron (MLP), Decision tree regressor
import argparse

# Refer to code ML_03_model_train.py and ML_07_boosting_models_test_run.ipynb

from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
import logging
import subprocess
import pandas as pd
import numpy as np
import os
import datetime

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
                        handlers=[logging.FileHandler(filename=log_fn, mode='w'), logging.StreamHandler()],
                        format='# %(name)s - %(levelname)s - %(message)s')

def parse_arguments():
    '''
    Parse arguments
    :return: Parsed arguments
    '''
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
    args = parser.parse_args()

    print('\n# Arguments used:')
    for arg in vars(args):
        print(f'# - {arg}:', getattr(args, arg))
    return parser



# #################### End of helper functions ####################

# Store some logging messages
msg = '#' * 80
msg += '\n# ' + '#'*20 + datetime.datetime.now().strftime('%Y-%m-%d') + '#'*20
msg += '\n#' * 80

args = parse_arguments()