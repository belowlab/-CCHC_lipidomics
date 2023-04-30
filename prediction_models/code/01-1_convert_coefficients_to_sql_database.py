# Convert coefficient files to SQL database
import sqlite3
import time
import argparse
import os

# #################### Helper functions ####################
def fill_coeff(lipid, coeffs, cur):
    '''
    Store coefficients of a given lipid into SQL database
    Params:
        - lipid: name of current lipid
        - coeffs: list of coefficients
        - cur: cursor of the SQL database
    Return:
        - count: Number of SNPs inserted (only include a SNP if coefficient is not 0)
    '''
    # Load GWAS stats file
    gwas_fn = f"{lipid.replace('(', '-').replace(')', '-').replace(' ', '_').replace('/', '-')}_SNPs_pval_0.001.txt"
    # print(gwas_fn, os.path.isfile(f'{args.gwas_result_dir}/{gwas_fn}'))
    fh = open(f'{args.gwas_result_dir}/{gwas_fn}')

    line = fh.readline()
    line = fh.readline().strip()
    count = 0
    while line != '':
        # GWAS summary stats has below format:
        # CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P
        # 1 chr1:828490:C:G	828490	G	C	2096	0.000715649	-2.02219	0.569577	0.000384725
        chr_num, snp_id, pos, _ = line.split(maxsplit=3)
        ref_allele, alt_allele = snp_id.split(':')[-2:]
        weight = coeffs.pop(0)
        if float(weight) != 0: # Only insert if coefficient is not zero
            # Table columns: chr, pos, snp_id, lipid, weight, ref_allele, eff_allele
            cur.execute(f"INSERT INTO weights VALUES ({chr_num}, {pos}, '{snp_id}', '{lipid}', {weight}, '{ref_allele}', '{alt_allele}')")
            count += 1
        line = fh.readline().strip()
    fh.close()
    return count

# #################### End of helper functions ####################

start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('--input',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/all_coeff_r2.txt',
                    help='input coefficient file')
parser.add_argument('--output_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params',
                    help='Output directory of SQL database')
parser.add_argument('--output_fn',
                    default='coeff.db',
                    help='Output file name of SQL database')
parser.add_argument('--gwas_result_dir',
                    default='/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3',
                    help='Directory of GWAS results')
args = parser.parse_args()
if args.output_dir.endswith('/'):
    args.output_dir = args.output_dir[:-1]
if args.gwas_result_dir.endswith('/'):
    args.gwas_result_dir = args.output_dir[:-1]

if os.path.isfile(f'{args.output_dir}/{args.output_fn}'):
    print(f'#File already exist: {args.output_dir}/{args.output_fn}')
    print('#Skip saving')
    exit()

print(f'#Create a SQl database at: {args.output_dir}/{args.output_fn}')
con = sqlite3.connect(f'{args.output_dir}/{args.output_fn}')
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
        lipid = line.split('\t')[0]
        coeffs = line.split()[-1].split(',')
        # Fill coefficients of current lipid into the database
        count = fill_coeff(lipid, coeffs, cur)
        print(f'# - Process {lipid}: {count} SNPs added')
        line = fh.readline().strip()
con.commit()
end_time = time.time()
print(f'#Finished in {(end_time-start_time)/60:.4f}min')