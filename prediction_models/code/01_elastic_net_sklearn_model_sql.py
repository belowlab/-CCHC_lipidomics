from sklearn.linear_model import ElasticNetCV
import sqlite3 # Try reading dosages from SQL database
import pandas as pd
import numpy as np
import os
import time
import argparse
import sys
# sys.path.append('/data100t1/home/wanying/lab_code/utils')
# from rank_based_inverse_normal_transformation import inverse_normal_transformation
import warnings
import datetime
warnings.filterwarnings(action='ignore')

'''
Example call:

# (Per AlexP) To run in base environment and avoid multi-threading conflict
# Run script as: OMP_NUM_THREADS=1 python my_script.py, for example:

# To run a single lipid
/usr/bin/time -v \
OMP_NUM_THREADS=1 python 01_elastic_net_sklearn_model_sql.py \
--output Sph\(d18:1\).txt \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params \
--dosage_db /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_train.db \
--multiallelic False \
--gwas_snps /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/\*_SNPs_pval_0.001.txt \
--lipid_name Sph\(d18:1\) \
--pheno /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted \
--n_alphas 100

# To run a list of lipids
OMP_NUM_THREADS=1 python 01_elastic_net_sklearn_model_sql.py \
--output lipid_species_1-10.txt.output \
--output_dir /data100t1/home/wanying/CCHC/lipidomics/prediction_models/elastic_net/training/model_params \
--dosage_db /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/sqlDB/dosage_train.db \
--multiallelic False \
--gwas_snps /data100t1/home/wanying/CCHC/lipidomics/output/traininig_set_lipid_species_GWAS/adj_for_sex_age_pval_1e-3/\*_SNPs_pval_0.001.txt \
--lipid_list /data100t1/home/wanying/CCHC/lipidomics/prediction_models/code/lipid_list/lipid_species_1-10.txt \
--pheno /data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted \
--n_alphas 100

'''
# ---------------------- Help functions ----------------------
def progress_bar(total, progress):
    '''
    Parameters:
    - total: total number of operations to be done
    - progress: number of operations completed
    '''
    percent = 100 * (progress/total)
    bar = '=' * int(percent) + '-' * int(100 - percent)
    print(f'|{bar}| {percent:.2f}%', end='\r')

def get_dosage_sql(dosage_cur, lst_snps):
    '''
        Get dosages from a SQL database based the provided lipid
        Param:
         - dosage_cur: cursor to the dosage database
         - lst_snps: a list of tuples of (chromosome number, SNP positions, REF allele, ALT allele) to be searched for
        Return:
         - sample_ids: IDs of genotype samples
         - dosage_matrix: dosage of given SNPs as a numpy array. Fill with NA if a SNP is not found
        '''
    # Index of column where dosage starts
    sample_id_indx = 6 # Columns in the database are: CHROM, POS, ID, REF, ALT, INFO, sample1, sample2, ...
    sample_ids = [x[0] for x in dosage_cur.execute('SELECT * FROM dosage LIMIT 1').description][sample_id_indx:]
    total_count, find, not_found = 0, 0, 0 # Track number of SNPs loaded and number of SNPs not found
    dosage_matrix, snps_loaded = [], [] # Store dosage and SNPs loaded (ie. found)
    print('# - Loading dosage:', len(lst_snps), 'SNPs')
    for snp in lst_snps:
        chr_num, pos, ref_allele, alt_allele = snp
        # print(f"\tSELECT * FROM dosage WHERE CHROM={chr_num} AND POS={pos} AND REF='{ref_allele}' AND ALT='{alt_allele}' LIMIT 1")
        result = dosage_cur.execute(f"SELECT * FROM dosage WHERE CHROM={chr_num} AND POS={pos} AND REF='{ref_allele}' AND ALT='{alt_allele}' LIMIT 1").fetchone()
        if result: # Append to output if not empty line
            dosage_matrix.append(result[sample_id_indx:])
            snps_loaded.append(snp)
            find += 1
        else:
            not_found += 1
        total_count += 1

        if total_count % 50 == 0:
            progress_bar(total=len(lst_snps), progress=total_count)
    progress_bar(total=len(lst_snps), progress=total_count)

    print(f'# - {total_count} SNPs checked, {find} SNPs loaded, {not_found} not found')
    # return sample_ids, np.array(dosage_matrix).reshape(-1, len(sample_ids))
    return sample_ids, np.array(dosage_matrix)

# Load dosage of all filtered SNPs by given p value threshold from GWAS
def load_all_dosage(gwas_snps, multiallelic, dosage_cur):
    '''
    Get dosage of SNPs from single-chrosmosome dosage files of a given lipid
    Params:
        - gwas_snps: Directory and file name of SNPs filtered by GWAS p value threshold
        - multiallelic: Drop multiallelic sites if False
        - dosage_cur: cursor to the SQL database of dosages
    Return:
        - start_time: start time of loading dosage
        - df_gwas_snp: a dataframe of GWAS SNPs
        - dosage_all: A numpy array of dosage. Each row is a SNP, each column is a subject
    '''
    # Check if file exists
    if not os.path.isfile(gwas_snps):
        print(f'# ERROR: GWAS SNP file not find: {gwas_snps}\n# END')
        exit()

    df_gwas_snp = pd.read_csv(gwas_snps, sep='\t').sort_values(by=['CHR', 'POS'])
    if not multiallelic: # Drop multiallelic SNPs
        df_gwas_snp.drop_duplicates(subset=['CHR', 'POS'], keep=False, inplace=True)
    df_gwas_snp['REF'] = df_gwas_snp['SNP'].apply(lambda x: x.split(':')[-2])
    df_gwas_snp['ALT'] = df_gwas_snp['SNP'].apply(lambda x: x.split(':')[-1])

    print('\n# Get dosage of GWAS SNPs to include in regression models')
    print('# - Checking by chromosome:')

    dosage_all = '' # A numpy array to store dosage from all chromosome
    start_time = time.time() # Time execution time
    # dosage_con = sqlite3.connect(dosage_db)
    # dosage_cur = dosage_con.cursor()
    for chr_num, df in df_gwas_snp.groupby(by='CHR'):
        print(f'# chr{chr_num}')
        lst_snps = list(df[['CHR', 'POS', 'REF', 'ALT']].itertuples(index=False, name=None)) # a list of tuples of (CHR, POS, REF, ALT)
        sample_ids, dosage_matrix = get_dosage_sql(dosage_cur=dosage_cur, lst_snps=lst_snps)
        if len(dosage_all) == 0: # if dosage array is empty (ie. the first chromosome)
            dosage_all = dosage_matrix
        else:
            dosage_all = np.append(dosage_all, dosage_matrix, axis=0)

    end_time = time.time()
    print(f'# - Checking finished in {(end_time-start_time):.4f}s')
    print('-' * 50)
    # dosage_con.close()
    return start_time, df_gwas_snp, dosage_all.astype('float64')
# ---------------------- End of help functions ----------------------

# ################# Process args #################
parser = argparse.ArgumentParser(description='Fit elastic net regression with 10 fold cross-validation')
parser.add_argument('-o', '--output', type=str,
                    help='Output file name to save alpha, l1_ratio and coefficients of chosen model')
parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='.')
parser.add_argument('--dosage_db', type=str, default='',
                    help='SQL database file of dosages')
parser.add_argument('--multiallelic', type=str,
                    default='False',
                    help='If false, multiallelic SNPs will be removed from model fitting')
parser.add_argument('--lipid_name', type=str, default='',
                    help='Name of a single lipid to be processed')
parser.add_argument('--lipid_list', type=str, default='',
                    help='A headerless file of lipids to be processed. Each line contains a single lipid. If provided, --lipid_name will be ignored')
parser.add_argument('--gwas_snps', type=str,
                    help='File name of SNPs filtered by GWAS p value threshold. Replace lipid name with * if --lipid_list. Eg. *_SNPs_pval_0.001.txt')
parser.add_argument('--pheno', type=str,
                    default='/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted',
                    help='Phenotype residuals adjusted by covariates such as age, sex, PCs')
parser.add_argument('--n_alphas', type=int, default=100,
                    help='Define how many alphas to test in CV. Dafault is 10. JTI used 100 as defined in R glmnet()')
args = parser.parse_args()
args.output = f"{args.output}.{datetime.datetime.now().strftime('%Y%m%d_%H:%M:%S')}"
if args.output_dir.endswith('/'): args.output_dir = args.output_dir[:-1]
if args.multiallelic.upper()[0]=='F' or args.multiallelic.upper()=='0':
    args.multiallelic = False
else: # Do not drop multiallelic sites if True
    args.multiallelic = True
print('# Run starts:', datetime.datetime.now().strftime('%Y-%m-%d'))
print('# Arguments used:')
for arg in vars(args):
    print(f'# - {arg}:', getattr(args, arg))

# Sanity checks
# Check if file exist
if not os.path.isfile(args.dosage_db):
    print(f'# ERROR: Dosage database not find: {dosage_db}\n# END')
    exit()
if '*' not in args.gwas_snps:
    print(f'# ERROR: --gwas_snps must replace lipid name with *\n# END')
    exit()

# ################# Load lipidomics data #################
print('# Load lipidomic data (lipid species)')
# fn_lipid = '/data100t1/home/wanying/CCHC/lipidomics/prediction_models/input_docs/lipid_traits_residuals/train/lipid_species_residuals_adj_for_sex_age_pc1-5.txt.reformatted'
df_lipid = pd.read_csv(args.pheno, sep='\t')
print(f"# - data loaded from {args.pheno.split('/')[-1]}: shape {df_lipid.shape}")

# Re-order lipidomic data so that sample IDs match the order in genotype file
fn_id_mapping = '/data100t1/home/wanying/CCHC/doc/samples_IDs/202211_merged_RNA_lipid_protein_genotype_mapping_and_availability.txt'
df_id_mapping = pd.read_csv(fn_id_mapping,
                            sep='\t').dropna(subset=['genotype_ID',
                                                     'lipidomic']).drop_duplicates(subset='lipidomic')[['LABID', 'genotype_ID']]

print(f'\n# Load genotype IDs for matching (only need to load the columns names of dosage database)')
dosage_con = sqlite3.connect(args.dosage_db)
dosage_cur = dosage_con.cursor()
res = dosage_cur.execute('SELECT * FROM dosage LIMIT 1')
df_genotype_id = pd.DataFrame([x[0] for x in res.description[6:]],
                              columns=['genotype_ID'])

print(f'# - Organize sample IDs so that their orders match in lipidomics data and dosage file')
df_lipid = df_genotype_id.merge(df_id_mapping.merge(df_lipid.drop_duplicates(subset='Sample ID'),
                                                    left_on='LABID',
                                                    right_on='Sample ID'), on='genotype_ID')
print(f'# - Final processed lipidomic data: {len(df_lipid)}')


# ################# Load GWAS snps of each lipid and run regression #################
# dosage_all: each row contains doages of a single SNP across all individuals
# !! Lip species PI(15-MHDA_20:4)\PI(17:0_20:4) is missing
# gwas_snp_dir = '/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_snps_pval_1e-3' # GWAS SNPs with p value<1e-3
# Save coefficients, alpha and l1 ratios of selected model for each lipid
output_fh = open(f'{args.output_dir}/{args.output}', 'w')
output_fh.write('lipid\talpha\tl1_ratio\tbest_r2\tcoefficients\n') # write header line
output_fh_lip_pred = open(f'{args.output_dir}/{args.output}.pred', 'w') # Save predicted values of each lipid using best fitted model
output_fh_lip_pred.write('Lipid'+'\t'+'\t'.join([val for val in df_lipid['Sample ID']])+'\n') # write header line

# Get list of lipids to be fitted
lst_lipid = []
if args.lipid_name != '':
    lst_lipid = [args.lipid_name]
elif args.lipid_list != '': # Create a list of lipids to be processed
    with open(args.lipid_list) as fh:
        line = fh.readline().strip()
        while line != '':
            lst_lipid.append(line)
            line = fh.readline().strip()

print(f'# {len(lst_lipid )} lipids to be processed')
for lipid_name in lst_lipid:
    count = 0
    # Modified lipid names
    lip = lipid_name.replace('(', '-').replace(')', '-').replace(' ', '_').replace('/', '-')
    gwas_snps = args.gwas_snps.replace('*', lip) # Fill * with modified lipid name
    if os.path.isfile(gwas_snps):
        # Get SNPs and dosage
        print(f'\n# Load GWAS SNPs for current lipid: {lipid_name}')
        load_dosage_start_time, df_gwas_snp, dosage_all = load_all_dosage(gwas_snps = gwas_snps,
                                                                          multiallelic = args.multiallelic,
                                                                          dosage_cur = dosage_cur)
        print(f'# - Number of SNPs loaded: {len(df_gwas_snp)}')

        print('# Run Elastic net regression')
        # lipid trait, already residuals and looks normal, so no need to INV
        y = df_lipid[lipid_name]
        # y = inverse_normal_transformation(df_lipid[lip])
        # print(y.shape)

        start_time = time.time()
        # Notes from sklearn docs:
        # - l1_ratio is the alpha in R glmnet
        # - alpha is the lambda in R gmlnet
        # Since PrediXcan used glmnet with alpha=0.5,and lambda selected by 10 fold cv,
        # The corresponding parameter in sklearn.ElasticNetCV() are:
        # - l1_ratio=0.5
        # - n_alphas=100, no user supllied selections for alpha, start with n_alphas=10 to save time (#TODO test how long it takes to run a full CV with 100 alphas)
        # - In R glmnet, when nobs > nvars, the default lambda.min.ratio is 0.0001
        # - 10 fold cv
        regr = ElasticNetCV(cv=10,
                            n_alphas=args.n_alphas,
                            random_state=0,
                            n_jobs=32,
                            l1_ratio=0.5) # Default l1 ratio=0.5
        X = dosage_all.T
        regr.fit(X, y)

        end_time = time.time()
        print(f'# - Model fitting finised in {(end_time - start_time):.4f}s')

        # Also output predicted values and best R2
        output_fh.write(f"{lipid_name}\t{regr.alpha_}\t{regr.l1_ratio_}\t{regr.score(X, y)}\t{','.join([str(x) for x in regr.coef_])}\n")
        output_fh_lip_pred.write(lipid_name+'\t'+'\t'.join([str(val) for val in regr.predict(X)])+'\n')
        print(f'# Total running time of current lipid: {(end_time - load_dosage_start_time)/60:.4f}m')
        # break
    else:
        print(f'# - Warning: {lipid_name} not found')
    count += 1
    print(f'# #################### {lipid_name}: {count} lipids processed processed ####################')
dosage_con.close()
output_fh.close()
output_fh_lip_pred.close()
