# Apply regression model on test set
# Modified from ~/CCHC/lipidomics/prediction_models/code/02_apply_model_with_tabix_NOT_DONE.py

# Apply model on test set (or other datasets)
'''
# Example call:
# Use this example!!!
lip_type=class
threshold=pval_1e-04_maf_0.01
python ML_06_apply_regression_model_on_test_set.py \
--vcf_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test/lipid_${lip_type} \
--vcf_files lipid_${lip_type}_chr*.pval_0.001_maf_0.05.test.vcf.gz \
--coeff_db /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/${lip_type}/merged_model_params/train_coeff_${lip_type}_${threshold}.db \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/predicted_values \
--output_fn_prefix test_${lip_type}_${threshold} \
--lipid_list /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}.list \
--overwrite False

# Dosage file function not work yet, use the vcf function to extract dosage
lip_type=class
threshold=pval_1e-04_maf_0.01
python ML_06_apply_regression_model_on_test_set.py \
--dosage_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/genotype_dosage/test/lipid_${lip_type} \
--dosage_files lipid_${lip_type}_chr*.pval_0.001_maf_0.05.test.vcf.dosage.gz \
--coeff_db /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/${lip_type}/merged_model_params/train_coeff_${lip_type}_${threshold}.db \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/predicted_values \
--output_fn_prefix test_${lip_type}_${threshold} \
--lipid_list /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_${lip_type}.list \
--overwrite False
'''

import pandas as pd
import sqlite3
import argparse
import datetime
import time
import numpy as np
import os
import logging
import subprocess

# ############## Helper functions ##############
def sanity_checks():
    '''
    Check input options
    :return:
    '''
    output_fn = os.path.join(args.output_dir, args.output_fn_prefix + '.pred')
    if os.path.isfile(output_fn) and not args.overwrite:
        logging.error('# - Output file already exists: ' + output_fn)
        logging.info('# - Skip saving. Exit')
        exit()
        
    # Check vcf or dosage files
    # --vcf_dir, --vcf_files, --dosage_dir, --dosage_files
    if (args.dosage_dir is not None) and (args.dosage_files is not None):
        logging.info('# Check dosage files')
        args.dosage_files = os.path.join(os.path.expanduser(args.dosage_dir), args.dosage_files)
        for i in range(1, 23): # Check chr1 to 22
            if not os.path.isfile(args.dosage_files.replace('*', str(i))):
                logging.error('# - Dosage file does not exist: ' + args.dosage_files.replace('*', str(i)))
                exit()
            if not os.path.isfile(args.dosage_files.replace('*', str(i))+'.tbi'):
                logging.error('# - Dosage index file does not exist: ' + args.dosage_files.replace('*', str(i)))
                exit()
        logging.info('# - PASS')
    elif (args.vcf_dir is not None) and (args.vcf_files is not None):
        logging.info('# Check vcf files')
        args.vcf_files = os.path.join(os.path.expanduser(args.vcf_dir), args.vcf_files)
        for i in range(1, 23):  # Check chr1 to 22
            if not os.path.isfile(args.vcf_files.replace('*', str(i))):
                logging.error('# - VCF file does not exist: ' + args.vcf_files.replace('*', str(i)))
                exit()
            if not os.path.isfile(args.vcf_files.replace('*', str(i))+'.tbi'):
                logging.error('# - VCF index file does not exist: ' + args.vcf_files.replace('*', str(i)))
                exit()
        logging.info('# - PASS')
    else:
        logging.error('# - No valid dosage file or vcf file was found; Exit')
        exit()

    # Check other input files
    if args.lipid_list is not None:
        logging.info('# Check lipid list')
        if not os.path.isfile(args.lipid_list):
            logging.error('# - Lipid list not found: ' + output_fn)
            logging.info('# Exit')
            exit()
        logging.info('# - PASS')
    elif args.lipid is not None:
        logging.info('# Search single lipid: ' + args.lipid)
        
def get_snp_weights(lipid, cur_coeff):
    '''
    Read from coefficient SQL database, load SNPS and weights of a given lipid
    Params:
        - lipid: name of the lipid
        - cur_coeff: cursor of coefficient database
    Return:
        - A dataframe with columns: 'chr', 'pos', 'snp_id', 'lipid', 'weight', 'ref_allele', 'alt_allele'.
                                    Weight is the column containing coefficients
    '''
    # The SQL database has a table "weights", with columns: 'chr', 'pos', 'snp_id', 'lipid', 'weight', 'ref_allele', 'alt_allele']
    res_coeff = cur_coeff.execute(f"SELECT * FROM weights WHERE lipid='{lipid}'")
    result_coeff = res_coeff.fetchall()
    df_coeff = pd.DataFrame(result_coeff, columns=['chr', 'pos', 'snp_id', 'lipid', 'weight', 'ref_allele', 'alt_allele'])

    # Return none if coefficient is not find
    if len(df_coeff) == 0:
        logging.info('# No weight was found in the SQL database for current lipid: '+lipid)
        return None
    else:
        return df_coeff


def get_dosage_from_vcf_single_chr(snp_pos, vcf_fn):
    '''
    Extract dosage of given SNPs from tabix indexed vcf file (single chromosome).
    Good for TOPmed imputed file. For vcfs from other imputation platform,
    may use get_dosage_from_processed_dosage_file() to avoid issues caused by file format.
    :param snp_pos: array-like of SNP positions, in format of chr1:123-123 (Use the same start and end for a single SNP)
    :param vcf_fn: Tabix indexed VCF file of single chromosome
    :return:
        - df_dosage: a dataframe of SNP dosages
        - c_found: number of SNPs found
        - c_multiallelic: number of multiallelic sites (droped)
    '''
    # Get header line. Only keep the last line as it is the column header
    # Get fields using cut: CHROM, POS, ID, REf, ALT and sample IDs
    cmd = f'tabix -H {vcf_fn} | cut -f 1-5,10-'
    res = subprocess.run(cmd, capture_output=True, shell=True, text=True).stdout.strip()
    col_header = res.split('\n')[-1].split()
    lst_all = [] # All values include CHROM, POS, ID, REf, ALT and dosage
    c_found, c_multiallelic = 0, 0 # Count total number of SNPs and multiallelic sites
    for pos in snp_pos:
        cmd = f'tabix {vcf_fn} {pos} | cut -f 1-5,10-'
        res = subprocess.run(cmd, capture_output=True, shell=True, text=True).stdout.strip()
        if len(res.split('\n'))>1: # Multiallelic sties
            # Drop multiallelic sties
            c_multiallelic += 1
            continue
        else:
            lst_dosage = [x.split(':')[1] for x in res.split('\t')[5:]]
            lst_other_cols = res.split('\t')[:5] # Non-dosage columns: CHROM, POS, ID, REf, ALT
            lst_all.append(lst_other_cols+lst_dosage)

        c_found += 1
    df_dosage = pd.DataFrame(data=lst_all, columns=col_header)
    return df_dosage, c_found, c_multiallelic

# Dosage file should be processed by ML_02_extract_snps_and_get_dosage.py and tabix indexed
# TODO: create a function to extract dosage from dosage files
def get_dosage_from_processed_dosage_file_single_chr(snp_pos, dosage_fn):
    pass

def calculate_pred_vals(lipid, cur_coeff):
    '''
    Calcualte predicted values of a single lipid
    :param lipid: Name of the lipid
    :param cur_coeff: cursor to coefficient database
    :return:
        - lst_sample_ids: list of sample IDs (for output file column header)
        - Predicted values of samples in genotype file
    '''
    df_coeff = get_snp_weights(lipid,
                               cur_coeff)  # Columns: 'chr', 'pos', 'snp_id', 'lipid', 'weight', 'ref_allele', 'alt_allele'
    if df_coeff is None:  # If lipid is not found, might due to naming format
        logging.info('# - Lipid not found, move to next\n')
        return None, None

    df_coeff['snp_pos'] = 'chr' + df_coeff['chr'].astype('str') + ':' + df_coeff['pos'].astype('str') + '-' + df_coeff[
        'pos'].astype('str')
    lst_df_dosage = []
    for chr_num, df in df_coeff.groupby(by='chr'):
        if args.dosage_files is not None:
            get_dosage_from_processed_dosage_file_single_chr()
        else:  # Use VCF files
            snp_pos, vcf_fn = df['snp_pos'], args.vcf_files.replace('*', str(chr_num))
            df_dosage, c_found, c_multiallelic = get_dosage_from_vcf_single_chr(snp_pos, vcf_fn)
            logging.info('# - Load dosage from CHR%s: Number of SNPs found/multiallelic sites ignored: %s/%s' % (
            chr_num, c_found, c_multiallelic))
        lst_df_dosage.append(df_dosage)
    df_dosage_all = pd.concat(lst_df_dosage)

    # Merge dosage with coefficient dataframe by chr number and pos
    # dosage file must have columns with column names:  '#CHROM', 'POS', 'ID', 'REF', 'ALT'
    # TODO: Check if ref and alt alleles are matched (necessary for external files, or files not imputed with TOPMed)
    assert ('POS' in df_dosage_all.columns) and ('#CHROM' in df_dosage_all.columns) and (
                'REF' in df_dosage_all.columns) and ('ALT' in df_dosage_all.columns)
    assert ('pos' in df_coeff.columns) and ('chr' in df_coeff.columns) and ('weight' in df_coeff.columns)
    # Get sample IDs
    lst_sample_ids = df_dosage_all.columns[5:]

    df_dosage_all['chr'] = pd.to_numeric(df_dosage_all['#CHROM'].apply(lambda x: int(x.split('chr')[-1])))
    df_dosage_all['POS'] = pd.to_numeric(df_dosage_all['POS'])
    df_coeff_dosage = df_coeff.merge(df_dosage_all, left_on=['chr', 'pos'], right_on=['chr', 'POS'])
    logging.info('# - Calculate predicted values')
    # Need to convert dosage values to float since they were read in as string
    pred_values = np.dot(df_coeff_dosage['weight'].values, df_coeff_dosage[lst_sample_ids].values.astype('float'))
    return lst_sample_ids, [str(val) for val in pred_values]  # Convert to string for output purpose

# ############## End of helper functions ##############

parser = argparse.ArgumentParser()
# Will check if dosage files are valid first then check vcf
# Use --dosage_dir, --dosage_files together
parser.add_argument('--dosage_dir',
                    default=None,
                    help='Directory to dosage file processed by ML_02_extract_snps_and_get_dosage.py')
parser.add_argument('--dosage_files',
                    default=None,
                    help='Dosage file processed by ML_02_extract_snps_and_get_dosage.py. Replace chr number with *. Will loop through 1 to 22')
# Use --vcf_dir, --vcf_files together.
parser.add_argument('--vcf_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_test_vcfs/',
                    help='Directory to genotype VCF files')
parser.add_argument('--vcf_files',
                    default='test_set_chr*.vcf.gz',
                    help='Genotype VCF file name. Replace chr number with *. Will loop through 1 to 22')

parser.add_argument('--coeff_db',
                    default='/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/prediction_models/elastic_net/class/merged_model_params/train_coeff_class_pval_1e-04_maf_0.01.db',
                    help='Coefficient SQL database file')
parser.add_argument('--output_dir',
                    default='./',
                    help='Output filename (including directory)')
parser.add_argument('--output_fn_prefix',
                    default='pred',
                    help='Prefix of output file')
parser.add_argument('--lipid',
                    default=None,
                    help='A single lipid to be processed.')
# lipid class list /data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/inputs/lipid_trait/lipid_class.list
parser.add_argument('--lipid_list',
                    default=None,
                    help='A list of lipid to test. One lipid per row without header line. Will ignore --lipid if this flag is provided and valid')
parser.add_argument('--overwrite', default='False', help='Overwrite output file if set to true')
args = parser.parse_args()

start_time = time.time()

# Convert overwrite to boolean value
if args.overwrite.upper() == 'FALSE' or args.overwrite == '0' or args.overwrite.upper().startswith('F'):
    args.overwrite = False
else: args.overwrite = True
    
# #################### Sanity checks ####################
args.output_dir = os.path.expanduser(args.output_dir)
msg = ''
if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)
    msg = f'# - Directory: not exist, create one at: {args.output_dir}'
    
log_fn = os.path.join(args.output_dir, args.output_fn_prefix+'.log')
print(f'# Log file saved at: {log_fn}')
logging.root.handlers = [] # Remove potential handler set up by others
# logging.info to console and save to log file. Change mode='a' for appending mode
logging.basicConfig(level=logging.INFO,
                    handlers=[logging.FileHandler(filename=log_fn, mode='w'), logging.StreamHandler()],
                    format='# %(name)s - %(levelname)s - %(message)s')
logging.info('Run start: ' + datetime.datetime.now().strftime('%Y-%m-%d') + '\n')


logging.info('# ' + '#' * 20 + ' Sanity checks ' + '#' * 20)
logging.info('# Check output directory and file')
if msg != '': logging.info(msg)
else: logging.info('# - PASS')

sanity_checks() # Check other inputs

logging.info('# Arguments passed:')
for arg, val in vars(args).items():
    logging.info(f'# - {arg}: {val}')
    
# #################### Load coefficients ####################
logging.info('\n')
logging.info('# '+ '#'*20 + ' Load coefficients ' + '#'*20)

# Connect to databases
logging.info('# Connect to coefficient databases')
con_coeff = sqlite3.connect(args.coeff_db)
cur_coeff = con_coeff.cursor()

logging.info('# Load lipid list')
fh_lipid_list = open(args.lipid_list)
lipid = fh_lipid_list.readline().strip()
count = 1
output_fn = os.path.join(args.output_dir, args.output_fn_prefix+'.pred')
fh_output = open(output_fn, 'w')
header_line = False # Track whether column names have been written to file
while lipid != '':
    msg = f'# Process lipid #{count}: {lipid}'
    logging.info(msg)
    sample_ids, pred_values = calculate_pred_vals(lipid, cur_coeff)

    if sample_ids is None: # If lipid not found, continue to next lipid
        count += 1
        lipid = fh_lipid_list.readline().strip()
        continue
        
    if not header_line: # Write header line if it has not been written yet
        fh_output.write('trait\t')
        fh_output.write('\t'.join(sample_ids) + '\n')
        header_line = True
    logging.info('# - Save to output file\n')
    fh_output.write(lipid + '\t')
    fh_output.write('\t'.join(pred_values)+'\n')
    lipid = fh_lipid_list.readline().strip()
    count += 1

fh_lipid_list.close()
con_coeff.close()
fh_output.close()

time_elapsed = time.time() - start_time
if time_elapsed > 60:
    logging.info('# Run finished in %0.4f min' % (time_elapsed/60))
else:
    logging.info('# Run finished in %0.4f sec' % (time_elapsed))

