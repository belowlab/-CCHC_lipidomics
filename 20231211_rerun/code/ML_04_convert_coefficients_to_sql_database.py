# Modified from ~/CCHC/lipidomics/prediction_models/code/01-1_convert_coefficients_to_sql_database.py
# Store regression coefficients to sql database for fast retrival

'''
Example call:
python /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/code/ML_04_convert_coefficients_to_sql_database.py \
--input /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/all_liplid_species_elastic_net_params.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/ \
--output_fn_prefix train_coeff \
--overwrite False

# Or below code if training outputs have not been merged into a single input file
threshold="pval_1e-07_maf_0.01"
python /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/code/ML_04_convert_coefficients_to_sql_database.py \
--search_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/class/${threshold} \
--model_param_file_suffix elastic_net \
--merged_model_param_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/class/${threshold} \
--merged_model_param_fn all_trait_merged_${threshold}.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/class/merged_model_params \
--output_fn_prefix train_coeff_${threshold} \
--overwrite False
'''

import sqlite3
import time
import argparse
import os
import logging

# #################### Helper functions ####################
def merge_model_parameter_files(input_dir, input_suffix, output_fn, output_dir=None):
    '''
    Search a folder to merge individual model parameter files into one
    :param input_dir: Directory to search individual files
    :param input_suffix: suffix of individual files (eg. .elastic_net)
    :param output_dir: Output directory. By default, it is the same as input_dir
    :param output_fn: Name of merged file
    :return: Save merged file using output_dir and output_fn. Return full path to saved file name, and number of files merged
    '''
    if output_dir is None: output_dir = input_dir
    output_fh = open(os.path.join(output_dir, output_fn), 'w')
    c = 0
    for fn in os.listdir(input_dir):
        if fn.endswith(input_suffix):
            with open(os.path.join(input_dir, fn)) as fh:
                # Assume the file has one header line and one result line
                # Fixed columns are: lipid (or other trait), alpha, l1_ratio, best_r2, coefficients, SNPs
                if c==0: output_fh.write(fh.readline())
                else: fh.readline()
                output_fh.write(fh.readline())
                c += 1
    logging.info('# - Merge individual parameter files into one: ' + str(c) + ' merged')
    output_fh.close()
    return os.path.join(output_dir, output_fn), c

def fill_coeff(lipid, coeffs, snps, cur):
    '''
    Store coefficients of a given lipid into SQL database
    Params:
        - lipid: name of current lipid
        - coeffs: list of coefficients
        - snps: list of snps used by the model
        - cur: cursor of the SQL database
    Return:
        - count: Number of SNPs inserted (only include a SNP if coefficient is not 0)
    '''
    count = 0
    for i in range(len(coeffs)):
        weight = coeffs[i]
        if float(weight) != 0: # Only insert if coefficient is not zero
            snp_id = snps[i] # SNP id in our dataset is chr:pos:ref:alt
            chr_num, pos, ref_allele, alt_allele = snp_id.split(':')
            chr_num = chr_num.split('chr')[-1]
            # Table columns: chr, pos, snp_id, lipid, weight, ref_allele, eff_allele
            cur.execute(f"INSERT INTO weights VALUES ({chr_num}, {pos}, '{snp_id}', '{lipid}', {weight}, '{ref_allele}', '{alt_allele}')")
            count += 1
    return count
# #################### End of helper functions ####################

start_time = time.time()
parser = argparse.ArgumentParser()
# Format of input coefficient file: row x column = trait x (lipid, alpha, l1_ratio, best_r2, coefficients, SNPs)
parser.add_argument('--input',
                    default='all_trait_elastic_net_params.txt',
                    help='Input coefficient file with columns named as: lipid (or other trait), alpha, l1_ratio, best_r2, coefficients, SNPs. Must has one header line and in this order')

# Use together: --search_dir, --model_param_file_suffix, --merged_model_param_dir and --merged_model_param_fn
# Once done, the merged file can be used as --input in future runs
parser.add_argument('--search_dir', default='',
                    help='Search directory of individual model parameter files and merge into one file. Will ignore --input if this flag is provided')
parser.add_argument('--model_param_file_suffix', default='', help='Suffix of individual model parameter files')
parser.add_argument('--merged_model_param_dir', default='', help='Directory of merged model parameter file. Will be tha same as --search_dir if leave empty')
parser.add_argument('--merged_model_param_fn', default='all_trait_elastic_net_params.txt', help='File name of merged model parameter file')

parser.add_argument('--output_dir',
                    default='./',
                    help='Output directory of SQL database')
parser.add_argument('--output_fn_prefix',
                    default='coeff',
                    help='Prefix of the output SQL database. ".db" will be added to the final file name.')
parser.add_argument('--overwrite',
                    default='False', choices=[0, 1, 'f', 't', 'F', 'T', 'False', 'True', 'false', 'true'],
                    help='If set to True, will overwrite existing output file with the same name')


args = parser.parse_args()

# Convert --overwrite to boolean value
if args.overwrite == '1' or args.overwrite.upper()[0] == 'T':
    args.overwrite = True
else: args.overwrite = False

args.input = os.path.expanduser(args.input)
args.search_dir = os.path.expanduser(args.search_dir)
args.merged_model_param_dir = os.path.expanduser(args.merged_model_param_dir)
args.output_dir = os.path.expanduser(args.output_dir)
output_fn = os.path.join(args.output_dir, args.output_fn_prefix+'.db')

log_fn = os.path.join(args.output_dir, args.output_fn_prefix+'.log')
print(f'\n# Log file saved at: {log_fn}')
logging.root.handlers = [] # Remove potential handler set up by others
# logging.info to console and save to log file. Change mode='a' for appending mode
logging.basicConfig(level=logging.INFO,
                    handlers=[logging.FileHandler(filename=log_fn, mode='w'), logging.StreamHandler()],
                    format='# %(name)s - %(levelname)s - %(message)s')

logging.info('#'*40+' Sanity checks '+'#'*40)
logging.info('# - Check input file or create one:')
if args.search_dir != '': # Create input file if --search_dir is provided
    logging.info('# - Search directory to create an input file')
    if not os.path.isdir(args.search_dir):
        logging.error('# - Search directory does not exist: '+args.search_dir + '\nEND')
        exit()
    else:
        if args.merged_model_param_dir == '':
            args.merged_model_param_dir = args.search_dir
        else:
            if not os.path.isdir(args.merged_model_param_dir):
                logging.info('# - Directory to save input file does not exist. Create one at ' + args.merged_model_param_dir )
                os.mkdir(args.merged_model_param_dir)
        # Create an input file and assign to --input flag
        args.input, c = merge_model_parameter_files(args.search_dir, args.model_param_file_suffix,
                                                    args.merged_model_param_fn, args.merged_model_param_dir)
        if c == 0: # If no individual model parameter file was merged, exit
            logging.error('# - No individual model parameter file was found. Exit')
            exit()
else: # Check if --input file exist
    if not os.path.isfile(args.input):
        logging.error('# - Input file does not exist. Create one with --model_param_file_suffix, --merged_model_param_dir and --merged_model_param_fn. Exit')
        exit()
    else: logging.info('# - PASS')

logging.info('# - Check output directory:')
if not os.path.isdir(args.output_dir):
    logging.info(f'# - Output directory does not exist, create one at: {args.output_dir}')
else: logging.info('# - PASS')

logging.info('# - Check output file: ')
if os.path.isfile(output_fn):
    logging.info(f'# - File already exist: {output_fn}')
    if not args.overwrite:
        logging.info('# - Skip saving.\nEND')
        exit()
    else: logging.info('# - Overwrite since --overwrite is True')
else: logging.info('# PASS')

logging.info('#'*40 + ' Arguments used ' + '#'*40)
for arg in vars(args):
    logging.info(f'# - {arg}: {getattr(args, arg)}')

logging.info('#'*40 + f' Create a SQl database at: {output_fn} ' + '#'*40)
# con = sqlite3.connect(f'{args.output_dir}/{args.output_fn}')
con = sqlite3.connect(output_fn)
cur = con.cursor()
cur.execute(f"DROP TABLE IF EXISTS weights") # Drop before create table
# Refer to JTI model columns: rsid, gene, weight, ref_allele, eff_allele
cur.execute("CREATE TABLE weights (chr INT, pos INT, snp_id TEXT, lipid TEXT, weight REAL, ref_allele TEXT, alt_allele TEXT)")

# Open coefficient file
logging.info('#'*40 + 'Load coefficient file' + '#'*40)
with open(args.input) as fh:
    line = fh.readline()
    line = fh.readline().strip()
    while line != '':
        lipid, alpha, l1_ratio, best_r2, coeffs, snps = line.split('\t')
        coeffs = coeffs.split(',')
        snps = snps.split(',')
        # Fill coefficientd and snpss used of current lipid into the database
        # Skip a SNP if coefficient is zero
        count = fill_coeff(lipid, coeffs, snps, cur)
        logging.info(f'# - Process {lipid}: {count} SNPs added')
        line = fh.readline().strip()
con.commit()
end_time = time.time()
logging.info('#'*40 + f' Finished in {(end_time-start_time)/60:.4f}min ' + '#'*40)