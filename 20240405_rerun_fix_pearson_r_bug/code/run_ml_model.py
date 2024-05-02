# Improved from code: /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/code/ML_03_model_train.py
'''
Example call: (remember to escape parenthesis in file name)

lip_type=species
lipid="PC(28:0)"
python run_ml_model.py \
--output_prefix "${lipid}" \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20240405_rerun_fix_pearson_r_bug/output/${lip_type} \
--train_vcf_path /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_${lip_type} \
--train_vcf_fn lipid_${lip_type}_chr*.pval_0.001_maf_0.05.vcf.gz \
--test_vcf_path /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test/lipid_${lip_type} \
--test_vcf_fn lipid_${lip_type}_chr*.pval_0.001_maf_0.05.test.vcf.gz \
--snp_fn /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_species_filter_by_pval_1e-05/PC-28:0-_SNPs.pval_1e-05.txt \
--trait_name ${lipid} \
--trait_fn /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}_ID_matched.no_dup.residual.all.txt \
--n_alphas 100 \
--drop_multiallelic True \
--reg_type elastic_net \
--id_col genotype_id \
--overwrite False
'''
from sklearn.metrics import mean_squared_error

from sklearn.linear_model import LassoCV, RidgeCV, ElasticNetCV
import pandas as pd
import numpy as np
import os
import argparse
import subprocess
import warnings
import time
import datetime
warnings.filterwarnings(action='ignore')
from scipy import stats
import logging
from joblib import dump # Save model


# Wanying's functions
from setup_logger import setup_log
from load_dosage import get_dosage_of_snps_from_file
from parse_args import parse_args

args = parse_args()

fn_model = os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.model.joblib')
if os.path.isfile(fn_model) and not args.overwrite:
    logging.info('# Output file exists. Set --overwrite True if really want to rerun')
    exit()
    
setup_log(os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.log'))
logging.info('#' * 80)
logging.info('# Run starts: %s' % datetime.datetime.now().strftime('%Y-%m-%d'))

logging.info('# Arguments used:')
for arg in vars(args):
    logging.info('# - %s: %s' % (arg, getattr(args, arg)))

# ################################### Load trait (lipidomics) data ##################################
msg = '# ' + '#'*20 + ' Load trait (lipidomics) data ' + '#'*20
logging.info(msg)
# Re-order lipidomics data so that sample IDs match the order in genotype file
# print(f'\n# Load genotype IDs for matching (only need to read the first line of dosage file)')
df_trait = pd.read_csv(args.trait_fn, sep='\t')
logging.info("# - Trait data loaded: shape %s*%s" % df_trait.shape)
logging.info('')
# print(df_trait.iloc[:5,:5])

# ################################## Load SNPs list of a given trait and run regression ##################################
# Get dosages
msg = '# ' + '#'*20 + ' Load dosages of SNPs of traning set: %s ' % args.trait_name + '#'*20
logging.info(msg)

start_time = time.time()
df_dosage = get_dosage_of_snps_from_file(args.train_vcf_path, args.train_vcf_fn, args.snp_fn, chr_col='CHR',
                                         pos_col='POS', id_col=args.id_col, drop_multiallelic=True)
end_time = time.time()
logging.info('# - Number of SNPs loaded: %s' % (df_dosage.shape[-1]-1)) # ID column does not count
logging.info('# - Variants lookup time: %.2f seconds' % (end_time - start_time))
logging.info('')
# print(df_dosage.iloc[:5, :5])

# Merge lipid trait with genotype data, prepare X and y for model training
df_merged = df_trait[[args.id_col, args.trait_name]].merge(df_dosage, on=args.id_col)
# print(df_merged.iloc[:5, :5])

y = df_merged[args.trait_name].values
X = df_merged.drop(columns=[args.id_col, args.trait_name])

logging.info('# ' + '#'*20 + ' Model training: %s regression' % args.reg_type + '#'*20)
logging.info('# - Number of samples into traininig: %s' % len(df_merged))
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

regr.fit(X, y)

end_time = time.time()
logging.info('# - Model training finished in %.2f seconds' % (end_time - start_time))
logging.info('')

#  ################################### Predict on testing set, calculate person r  ###################################
logging.info('# ' + '#'*20 + ' Load dosages of SNPs of test set ' + '#'*20)
df_dosage_test = get_dosage_of_snps_from_file(args.test_vcf_path, args.test_vcf_fn, args.snp_fn, chr_col='CHR',
                                              pos_col='POS', id_col=args.id_col, drop_multiallelic=True)
# Merge lipid trait with genotype data, prepare X and y for model training
df_merged_test = df_trait[[args.id_col, args.trait_name]].merge(df_dosage_test, on=args.id_col)
y_test = df_merged_test[args.trait_name].values
X_test = df_merged_test.drop(columns=[args.id_col, args.trait_name])
logging.info('# - Number of samples into testing: %s' % len(df_merged_test))
logging.info('')

# Calcualte pearon r on training and testing set
logging.info('# ' + '#'*20 + ' Model evaluation ' + '#'*20)
pred_train = regr.predict(X)
pred_test = regr.predict(X_test)

logging.info('# - Training set: pearson_r=%s, pval=%s' % stats.pearsonr(y, pred_train))
logging.info('# - Test set: pearson_r=%s, pval=%s' % stats.pearsonr(y_test, pred_test))



# ---------------------------- Log additional metrics ----------------------------
logging.info('# - Training set: spearman_r=%s, pval=%s' % stats.spearmanr(y, pred_train))
logging.info('# - Test set: spearman_r=%s, pval=%s' % stats.spearmanr(y_test, pred_test))

logging.info('# - Training set: mse=%s' % mean_squared_error(y, pred_train))
logging.info('# - Test set: mse=%s' % mean_squared_error(y_test, pred_test))

# --------------------------------------------------------------------------------




#  ################################### Save fitted model, coefficients, intercept and predicted values ###################################
# Save fitted model
model_fn = os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.model.joblib')
dump(regr, model_fn)

# Save coefficients and intercept
coefficients = np.insert(regr.coef_, 0, regr.intercept_) # Include intercept in the output
df_params = pd.DataFrame(data={'feature':np.insert(regr.feature_names_in_, 0, 'intercept'), 'weight':coefficients})
df_params.to_csv(os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.params'), sep='\t', index=False)

# Save predicted values (train and test)
df_pred_train = pd.DataFrame(data={args.id_col:df_merged[args.id_col], args.trait_name:pred_train})
df_pred_test = pd.DataFrame(data={args.id_col:df_merged_test[args.id_col], args.trait_name:pred_test})
df_pred_train.to_csv(os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.train.pred'), sep='\t', index=False)
df_pred_test.to_csv(os.path.join(args.output_dir, args.output_prefix+'.'+args.reg_type+'.test.pred'), sep='\t', index=False)

