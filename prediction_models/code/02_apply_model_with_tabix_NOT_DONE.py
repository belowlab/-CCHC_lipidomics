# Apply model on test set (or other datasets)
'''
# Example call:
python 02_apply_model.py \
--dosage /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db \
--coeff /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/train_coeff.db \
--output /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/test/test_pred_vals.txt \
--lipid_list /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_list.txt \
--overwrite True
'''

import pandas as pd
import sqlite3
import argparse
import datetime
import time
import numpy as np
import os

# ############## Helper functions ##############
def get_pred_vals(lipid, cur_coeff, cur_dosage, sample_ids):
    '''
    Read from SQL databases, load weights and dosage of SNPs of a given lipid. Calculate predicted values
    Params:
        - lipid: numbe of the lipid
        - cur_coeff: cursor of coefficient database
        - cur_dosage: cursor of dosage database
        - sample_ids: sample IDs
    Return:
        - n_missing: number of SNPs not found in dosage file
        - pred_vals: a numpy array of predicted values for the given lipid
    '''
    res_coeff = cur_coeff.execute(f"SELECT snp_id, weight FROM weights WHERE lipid='{lipid}'")
    result_coeff = res_coeff.fetchall()
    snp = ["'"+val[0]+"'" for val in result_coeff] # Get SNP IDs
    # coeff_array = np.array([val[1] for val in result_coeff]) # Load coefficients
    df_coeff = pd.DataFrame(result_coeff).rename(columns={0:'snp_id', 1:'weight'})
    res_dosage = cur_dosage.execute(f"SELECT {'ID,' + ','.join(sample_ids)} FROM dosage WHERE ID IN ({','.join(snp)})")
    df_dosage = pd.DataFrame(res_dosage.fetchall()).rename(columns={0:'snp_id'})
    
    # Return none if coefficient or dosage is not find
    if len(df_coeff)==0 or len(df_dosage)==0:
        return None, None
    
    n_missing = f'{len(df_coeff)-len(df_dosage)}/{len(df_coeff)}' # Number of SNPs missing in dosage file, in format as missing/total
    
    df_merged = df_coeff.merge(df_dosage, on='snp_id')
    
    # Calculate dot product of weight and dosage matrices
    pred_vals = np.dot(df_merged['weight'].values, df_merged.iloc[:, 2:].values)
    return n_missing, pred_vals
# ############## End of helper functions ##############


print('Last run:', datetime.datetime.now().strftime('%Y-%m-%d'))
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--dosage', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_test.db',
                    help='Input SQL database containing SNP dosage')
parser.add_argument('--coeff', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params/coeff_100_alpha.db',
                    help='Coefficient SQL database file')
parser.add_argument('--output', default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/test/pred_vals.txt',
                    help='Output filename (including directory)')
parser.add_argument('--lipid_list',
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_list.txt',
                    help='A list of lipid to test. One lipid per row without header line.')
parser.add_argument('--overwrite', default=False, help='When set to True, overwrite output file if it already exists')
args = parser.parse_args()

# print(args, type(args))
print('#Arguments passed:')
for arg, val in vars(args).items():
    print(f'# - {arg}: {val}')

if args.overwrite.upper()=='FALSE' or args.overwrite=='0' or args.overwrite.upper().startswith('F'): 
    args.overwrite = False
else:
    args.overwrite = True

if os.path.isfile(args.output) and not args.overwrite:
    print(f'#Output file already exists: {args.output}')
    print('#Skip saving. Exit')
    exit()
    
fh_output = open(args.output, 'w')

# Connect to databases
print('#Connect to dosage and coefficient databases')
con_coeff = sqlite3.connect(args.coeff)
cur_coeff = con_coeff.cursor()

con_dosage = sqlite3.connect(args.dosage)
cur_dosage = con_dosage.cursor()
# Get sample IDs
headers = cur_dosage.execute('SELECT * FROM dosage').description
sample_ids = [val[0] for val in headers][6:]

# Write header row to output file
fh_output.write('Lipid\t'+'\t'.join(sample_ids)+'\n')

print('#Load lipid list')
fh_lipid_list = open(args.lipid_list)
lipid = fh_lipid_list.readline().strip()

count = 1
while lipid != '':
    print(f'# - Process lipid #{count}: {lipid}', end='')
    n_missing, pred_vals = get_pred_vals(lipid, cur_coeff, cur_dosage, sample_ids)
    if pred_vals is None: # If the lipid is not found
        print('; Lipid not found')
    else:
        fh_output.write(f'{lipid}\t'+'\t'.join(pred_vals.astype('U'))+'\n')
        print(f'; {n_missing} SNPs not found in dosage database')
    
    lipid = fh_lipid_list.readline().strip()
    count += 1

fh_lipid_list.close()
con_coeff.close()
con_dosage.close()
fh_output.close()

end_time = time.time()
print(f'#Run finished in {(end_time-start_time)/60:.4f} min')

