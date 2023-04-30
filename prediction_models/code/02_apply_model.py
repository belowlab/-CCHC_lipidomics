import sqlite3
import argparse
import datetime
import time
import numpy as np

# ############## Helper functions ##############
def get_dosage(lipid):
    '''
    Read from SQL database and load dosage of a given lipid
    Return: 1D numpy array of dosage
    '''
    gwas_fn = f'{lipid}_SNPs_pval_0.001.txt' # File name of GWAS summary stats
    pass
# ############## End of helper functions ##############


print('Last run:', datetime.datetime.now().strftime('%Y-%m-%d'))
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--dosage', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db',
                    help='Input SQL database containing SNP dosage')
parser.add_argument('--coeff', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/all_coeff_r2.txt',
                    help='Coefficient file')
parser.add_argument('--gwas_dir', default='/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3',
                    help='Coefficient file')
args = parser.parse_args()

if args.gwas_dir.endswith('/'):
    args.gwas_dir = args.gwas_dir[:-1]

# Load coefficients
with open(args.coeff) as fh:
    line = fh.readline()
    line = fh.readline().strip()
    while line != '':
        lipid = line.split()[0]
        coeffs = line.split()[-1].split(',')
        print(f'#Process lipid: {lipid}')
        # Count number of variants used in the model
        num_snps_in_model = 0
        coeffs_lst = []
        for val in coeffs:
            if float(val) != 0: num_snps_in_model += 1
            coeffs_lst.append(float(val))
        # print(len(coeffs), ':', num_snps_in_model)
        line = fh.readline()
        coeffs_array = np.array(coeffs_lst) # Coefficients of current lipid in numpy array
        # Load SNP dosage
        dosage = get_dosage(lipid)