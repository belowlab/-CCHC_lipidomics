import pandas as pd
import numpy as np
import os
import subprocess

try:
    from tqdm import tqdm
except:
    print('# install package tqdm')
    subprocess.run('pip install tqdm'.split())

def get_dosage(vcf_fn, pos, get_sample_ids=False):
    '''
    Get dosages of a list of given positions
    Params
    - vcf_fn: file name of the vcf to search
    - pos: an array like positions to lookup. For example: ['chr22:10513905-10513905', 'chr22:10514201-10514201']
    - get_sample_ids: will get sample ids and return if True
    Return
    - sample_ids: a list sample IDs from vcf
    - snp_ids: ids of SNPs with dosage found
    - dosages: a matrix of dosages in ist form. for example: [[0,1,1], [0,0,0]]
    '''
    chromosome = pos[0].split(':')[0] # Get chromosome number for loggging purpose
    sample_ids, snp_ids = [], []
    
    if get_sample_ids:
        # Get the header line containing sample IDs
        cmd_header = f'tabix -H {os.path.join(vcf_path, vcf_fn)}' 
        sample_ids = subprocess.run(cmd_header.split(),
                                    capture_output=True,
                                    text=True).stdout.rstrip().split('\n')[-1].split(maxsplit=9)[-1].split()
    # Get dosages of all SNPs
    pos_str = ' '.join(pos)
    cmd = f'tabix {os.path.join(vcf_path, vcf_fn)} {pos_str}' # Find all SNPs at once
    lines = subprocess.run(cmd.split(), capture_output=True, text=True).stdout
    dosages = []
    if len(lines) != 0: # If found something, lines will not be empty
        for line in tqdm(lines.split('\n'), desc=f'# - Getting dosage from {chromosome}'):
            try: # Last item will be a empty string, use try-except to avoid
                dosage = [float(val.split(':')[1]) for val in line.split(maxsplit=9)[-1].split()]
                dosages.append(dosage)
                snp_ids.append(line.split(maxsplit=3)[2])
            except:
                continue
    return sample_ids, snp_ids, dosages

def get_dosage_of_snps_from_gwas_output(vcf_path, vcf_fn, snp_fn):
    '''
    Given a (filtered) GWAS output file, get dosages of all SNPs from the file
    Params
    - vcf_path, vcf_fn: path and file name to the vcfs. vcf_fn must have chromosome number replaced with *, such as 'max_unrelated_set_chr*.vcf.gz'
    - snp_fn: file name to the filtered SNPs, such as '/data100t1/home/wanying/CCHC/lipidomics/TG-123.txt'
    Return
    Return None if CHR and POS are not columns in the SNP file.
    Otherwise, return a Dataframe of dosage with sample IDs and SNP IDs labeled.
    Results have shape of sample x features(SNPs)
    
    # Path and file name examples
    vcf_path='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs'
    vcf_fn = 'max_unrelated_set_chr*.vcf.gz' # Replace * with chromosome number in the loop
    snp_path = '/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_species_filter_by_pval_1e-05'
    snp_fn = 'TG-58:8-_[NL-22:6]_SNPs.pval_1e-05.txt'
    '''
    
    df_snps = pd.read_csv(snp_fn, sep='\t', dtype='str')
    
    # 'CHR' and 'POS' columns shoulbe found in df_snps.columns
    if not 'CHR' in df_snps.columns:
        print('# Error: Column CHR not found in the filtered SNP file')
        return
    if not 'POS' in df_snps.columns:
        print('# Error: Column CHR not found in the filtered SNP file')
        return
    
    # Get list of positions from the SNP dataframe
    df_snps['position_for_tabix'] = 'chr' + df_snps['CHR'] + ':' + df_snps['POS'] + '-' + df_snps['POS']
    dosage_all, snp_ids_all = [], [] # Store dosages from each chromosome

    for i, (chr_num, df) in enumerate(df_snps.groupby('CHR')):
        positions = df['position_for_tabix'].values
        if i==0:
            sample_ids, snp_ids, dosages = get_dosage_matrix(vcf_fn=os.path.join(vcf_path, vcf_fn.replace('*', chr_num)),
                                                             pos=positions,
                                                             get_sample_ids=True)
        else:
            _, snp_ids, dosages = get_dosage_matrix(vcf_fn=os.path.join(vcf_path, vcf_fn.replace('*', chr_num)),
                                                             pos=positions)
        dosage_all += dosages
        snp_ids_all += snp_ids

    # df_dosage = pd.DataFrame(data=np.array(dosage_all), columns=header, index=snp_ids_all).T
    return pd.DataFrame(data=np.array(dosage_all), columns=sample_ids, index=snp_ids_all).T