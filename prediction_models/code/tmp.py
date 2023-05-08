import sqlite3
import argparse
import datetime
import time
import numpy as np
import os

print('Last run:', datetime.datetime.now().strftime('%Y-%m-%d'))
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--dosage', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db',
                            help='Input SQL database containing SNP dosage')
parser.add_argument('--coeff', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/coeff_100_alpha.db',
                            help='Coefficient SQL database file')
parser.add_argument('--output', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/test/pred_vals.txt',
                            help='Output filename (including directory)')
parser.add_argument('--lipid_list', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_list.txt',
                            help='A list of lipid to test. One lipid per row without header line.')
args = parser.parse_args('--dosage /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db'.split())

# if os.path.isfile(args.output):
#     print(f'#Output file already exists: {args.output}')
#     print('#Skip saving. Exit')
#     exit()

# Connect to databases
print('#Connect to dosage and coefficient databases')
con_coeff = sqlite3.connect(args.coeff)
cur_coeff = con_coeff.cursor()

con_dosage = sqlite3.connect(args.dosage)
cur_dosage = con_dosage.cursor()
# Get sample IDs
headers = cur_dosage.execute('SELECT * FROM dosage').description
sample_ids = [val[0] for val in headers][6:]

print('#Load lipid list')
fh_lipid_list = open(args.lipid_list)
lipid = fh_lipid_list.readline().strip()
count = 1
# while lipid != '':
#     print(f'# - Process lipid #{count}: {lipid}')

#     pred_vals = get_pred_vals(lipid, cur_coeff, cur_dosage, sample_ids)
#     lipid = fh_lipid_list.readline().strip()
#     count += 1

res_coeff = cur_coeff.execute(f"SELECT snp_id, weight FROM weights WHERE lipid='{lipid}'")
result_coeff = res_coeff.fetchall()
snp = ["'"+val[0]+"'" for val in result_coeff] # Load SNP IDs
coeff_array = np.array([val[1] for val in result_coeff]) # Load coefficients

# Order is retained as original order in the database. So SNPs and coefficients are lined up already, no need to sort
res_dosage = cur_dosage.execute(f"SELECT {','.join(sample_ids)} FROM dosage WHERE ID IN ({','.join(snp)})")
dosage = res_dosage.fetchall()

end_time = time.time()
print(f'#Run finished in {(end_time-start_time)/60:.4f} min')


# fh_lipid_list.close()
# con_coeff.close()
# con_dosage.close()print('Last run:', datetime.datetime.now().strftime('%Y-%m-%d'))
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--dosage', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db',
                            help='Input SQL database containing SNP dosage')
parser.add_argument('--coeff', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/coeff_100_alpha.db',
                            help='Coefficient SQL database file')
parser.add_argument('--output', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/test/pred_vals.txt',
                            help='Output filename (including directory)')
parser.add_argument('--lipid_list', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_list.txt',
                            help='A list of lipid to test. One lipid per row without header line.')
args = parser.parse_args('--dosage /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db'.split())

# if os.path.isfile(args.output):
#     print(f'#Output file already exists: {args.output}')
#     print('#Skip saving. Exit')
#     exit()

# Connect to databases
print('#Connect to dosage and coefficient databases')
con_coeff = sqlite3.connect(args.coeff)
cur_coeff = con_coeff.cursor()

con_dosage = sqlite3.connect(args.dosage)
cur_dosage = con_dosage.cursor()
# Get sample IDs
headers = cur_dosage.execute('SELECT * FROM dosage').description
sample_ids = [val[0] for val in headers][6:]

print('#Load lipid list')
fh_lipid_list = open(args.lipid_list)
lipid = fh_lipid_list.readline().strip()
count = 1
# while lipid != '':
#     print(f'# - Process lipid #{count}: {lipid}')

#     pred_vals = get_pred_vals(lipid, cur_coeff, cur_dosage, sample_ids)
#     lipid = fh_lipid_list.readline().strip()
#     count += 1

res_coeff = cur_coeff.execute(f"SELECT snp_id, weight FROM weights WHERE lipid='{lipid}'")
result_coeff = res_coeff.fetchall()
snp = ["'"+val[0]+"'" for val in result_coeff] # Load SNP IDs
coeff_array = np.array([val[1] for val in result_coeff]) # Load coefficients

# Order is retained as original order in the database. So SNPs and coefficients are lined up already, no need to sort
res_dosage = cur_dosage.execute(f"SELECT {','.join(sample_ids)} FROM dosage WHERE ID IN ({','.join(snp)})")
dosage = res_dosage.fetchall()

end_time = time.time()
print(f'#Run finished in {(end_time-start_time)/60:.4f} min')


# fh_lipid_list.close()
# con_coeff.close()
# con_dosage.close()


missing = []
for s in snp:
    number_find = len(cur_dosage.execute(f"SELECT ID FROM dosage WHERE ID={s}").fetchall())
    if number_find==0:
        missing.append(f'{s}|{number_find}')
    else: print('.', end='', flush=True)
print('################# Missing #################')
for val in missing:
    print(val)
