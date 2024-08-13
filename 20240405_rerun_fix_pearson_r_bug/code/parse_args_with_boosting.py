# Process arguments and sanity checks
import os
import time
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Fit regression model of choice (elastic net, ridge, lasso or OLS)')
    parser.add_argument('-o', '--output_prefix', type=str,
                               help='Output file to save alpha, l1_ratio and coefficients of chosen model')
    parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='.')
    parser.add_argument('--train_vcf_path', type=str, help='Derictory to vcf files for model training',
                        default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/train/lipid_species')
    parser.add_argument('--train_vcf_fn', type=str, help='File name format of the vcf files for model training. Use * to replace chromosome number',
                        default='species_chr*.vcf.gz')
    parser.add_argument('--test_vcf_path', type=str,
                        help='Derictory to vcf files for testing. Use together with --test_vcf_fn',
                        default=None)
    parser.add_argument('--test_vcf_fn', type=str,
                        help='File name format of the vcf files for testing. Use * to replace chromosome number. Will test performance if provided',
                        default=None)
    parser.add_argument('--snp_fn', type=str,
                        help='File name to a list of SNPs (eg. GWAs SNPs with pval<1e-3). Must contains chromosome and position columns. Must be a tsv file',
                        default='AC-10:0-_SNPs_pval_0.001.txt')
    parser.add_argument('--trait_name', type=str,
                        help='Name of the trait (lipid) in the --trait_fn to be processed.')
    parser.add_argument('--trait_fn', type=str,
                        help='File name of lipidomics data (or other trait). File must in sample x trait format. Must have a column matches genotype IDs. Contain both train and test data if need evaluation on the test set',
                        default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_species_ID_matched.no_dup.residual.train.txt')
    parser.add_argument('--n_alphas', type=int, default=100,
                        help='For elastic net. Defines how many alphas to test in CV. JTI used 100 as defined in R glmnet()')
    parser.add_argument('--drop_multiallelic', type=str, default='True',
                        help='If true, multiallelic SNPs will be removed from model fitting')
    parser.add_argument('--model', type=str, default='elastic_net', choices=['elastic_net', 'ridge', 'lasso', 'ols', 'Ada', 'Gradient'],
                        help="Type of regression. Choose from: 'elastic_net', 'ridge', 'lasso' and 'ols'")
    parser.add_argument('--id_col', type=str, default='genotype_id',
                        help="Name of the ID column in trait file. Default is genotype_id")
    parser.add_argument('--overwrite', type=str, default='False',
                        help='If true, will rerun and overwrite existing files')
    args = parser.parse_args()

    # Sanity checks
    sanity_checks(args)
    
    return args

def sanity_checks(args):
    # Make all directory to absolute path
    args.output_dir = os.path.expanduser(args.output_dir)
    args.trait_fn = os.path.expanduser(args.trait_fn)

    # Check if files exist
    print('# - Check SNP file: ', end='') # (filtered) GWAS SNP
    if not os.path.isfile(args.snp_fn):
        print('A list of SNPs not found:', args.snp_fn)
        exit()
    else:
        print('PASS')

    print('# - Check trait file (eg. residuals of lipidomic measures): ', end='')
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

    if args.drop_multiallelic.upper()[0]=='F' or args.drop_multiallelic=='0':
        args.drop_multiallelic = False
    else: # Do not drop multiallelic sites if True
        args.drop_multiallelic = True
        
    if args.overwrite.upper()[0]=='F' or args.overwrite=='0':
        args.overwrite = False
    else:
        args.overwrite = True