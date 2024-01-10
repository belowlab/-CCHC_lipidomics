# Convert coefficient files to SQL database
'''
Example call:
python 01-1_convert_coefficients_to_sql_database.py \
--input /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/all_liplid_species_elastic_net_params.txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/ \
--output_fn train_coeff.db
'''

import sqlite3
import time
import argparse
import os

# #################### Helper functions ####################
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
parser.add_argument('--input',
                    default='all_liplid_species_elastic_net_params.txt',
                    help='input coefficient file with columns: lipid, alpha, l1_ratio, best_r2, coefficients, SNPs')
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/',
                    help='Output directory of SQL database')
parser.add_argument('--output_fn',
                    default='coeff.db',
                    help='Output file name of SQL database')

args = parser.parse_args()
if args.output_dir.endswith('/'):
    args.output_dir = args.output_dir[:-1]

if os.path.isfile(f'{args.output_dir}/{args.output_fn}'):
    print(f'#File already exist: {args.output_dir}/{args.output_fn}')
    print('#Skip saving')
    exit()

print('#Arguments used:')
for arg in vars(args):
    print(f'# - {arg}: {getattr(args, arg)}')    

print(f'\n#Create a SQl database at: {args.output_dir}/{args.output_fn}')
# con = sqlite3.connect(f'{args.output_dir}/{args.output_fn}')
con = sqlite3.connect(os.path.join(args.output_dir, args.output_fn))
cur = con.cursor()
cur.execute(f"DROP TABLE IF EXISTS weights") # Drop before create table
# Refer to JTI model columns: rsid, gene, weight, ref_allele, eff_allele
cur.execute("CREATE TABLE weights (chr INT, pos INT, snp_id TEXT, lipid TEXT, weight REAL, ref_allele TEXT, alt_allele TEXT)")

# Open coefficient file
print('#Load coefficient file')
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
        print(f'# - Process {lipid}: {count} SNPs added')
        line = fh.readline().strip()
con.commit()
end_time = time.time()
print(f'#Finished in {(end_time-start_time)/60:.4f}min')