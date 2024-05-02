import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import argparse
import subprocess
import warnings
import logging
warnings.filterwarnings(action='ignore')

def get_dosage_matrix(vcf_fn, pos, get_sample_ids=False, indx_dosage_field=1, drop_multiallelic=True):
    '''
    Get dosages of a list of given positions from a bgzipped and tabixed genotype vcf file (single chromosome VCF)
    Params
    - vcf_fn: file name of the vcf to search
    - pos: an array like positions to lookup. For example: ['chr22:10513905-10513905', 'chr22:10514201-10514201']
    - get_sample_ids: will get sample ids and return if True
    - indx_dosage_field = 1: Index of dosage field. For GT:DS:HDS:GP format the value is 1 (eg. in TOPmed imputed vcf)
    Return
    - header: a list sample IDs from vcf
    - dosages: a matrix of dosages in ist form. for example: [[0,1,1], [0,0,0]]
    - drop_multiallelic: ignore multiallelic sites if true
    '''
    chromosome = pos[0].split(':')[0] # Get chromosome number for logging purpose
    sample_ids, snp_ids, dosages = [], [], []
    
    if get_sample_ids:
        # Get the header line containing sample IDs
        cmd_header = f'tabix -H {vcf_fn}' 
        sample_ids = subprocess.run(cmd_header.split(),
                                    capture_output=True,
                                    text=True, check=True).stdout.rstrip().split('\n')[-1].split(maxsplit=9)[-1].split()
    
    if not drop_multiallelic:
        # Get dosages of all SNPs at once (tabix allows chaining of multiple positions)
        pos_str = ' '.join(pos)
        cmd = f'tabix {vcf_fn} {pos_str}' # Find all SNPs at once
        lines = subprocess.run(cmd.split(), capture_output=True, text=True, check=True).stdout.rstrip() # Remove last newline character
        if len(lines) != 0: # If found something, lines will not be empty
            for line in tqdm(lines.split('\n'), desc=f'# - Getting dosage from {chromosome}'):
                # TOPmed imputed VCF has format GT:DS:HDS:GP. Will not work if dosage (DS) field is not the second one (need to manually fix)
                dosage = [float(val.split(':')[indx_dosage_field]) for val in line.split(maxsplit=9)[-1].split()]
                dosages.append(dosage)
                snp_ids.append(line.split(maxsplit=3)[2])
        msg = '# CHR%s: keep multiallelic sites' % chromosome
    else:
        # Lookup SNPs one by one. Skip a SNP if multiple lines returned
        count_multiallelic, count_not_found = 0, 0 # Record number of multiallelic sites and missing SNPs
        for position in tqdm(pos, desc=f'# - Getting dosage from {chromosome}'):
            cmd = f'tabix {vcf_fn} {position}' # Find one SNPs at a time
            line = subprocess.run(cmd.split(), capture_output=True, text=True, check=True).stdout.rstrip() # Remove last newline character
            if len(line.split('\n')) == 1:
                if line == '': # SNP not found
                    count_not_found += 1
                    continue
                dosage = [float(val.split(':')[indx_dosage_field]) for val in line.split(maxsplit=9)[-1].split()]
                dosages.append(dosage)
                snp_ids.append(line.split(maxsplit=3)[2])
            elif len(line.split('\n')) > 1: # Multiallelic site
                count_multiallelic += 1
                continue
        msg = '# CHR%s: Number of multiallelic sites dropped=%s; Number of SNPs not found=%s' % (chromosome, count_multiallelic , count_not_found)
    logging.info(msg)
    return sample_ids, snp_ids, dosages

def get_dosage_of_snps_from_file(vcf_path, vcf_fn, snp_fn, chr_col='CHR', pos_col='POS', id_col='genotype_id', drop_multiallelic=True):
    '''
    Given a (filtered) GWAS output file, get dosages of all SNPs from the VCF files
    Params:
    - vcf_path, vcf_fn: path and file name to the vcfs. vcf_fn must have chromosome number replaced with *, such as 'max_unrelated_set_chr*.vcf.gz'
    - snp_fn: name of the SNP list file. The format is a tab delimited (tsv) file with one SNP per line. Chromosome and position columns must be numeric values.
              Chromosome and position is indicated by the "CHR" and "POS" columns. Other columns will not be used.
              Example: (such as fastGWA output)
              CHR  SNP                  POS    A1  A2      N     AF1          BETA        SE        P
                1    chr1:1234:G:A       1234  A   G       1606  0.002   0.1    0.3  0.5
                2    chr2:1234:G:A       1234  A   G       1606  0.0007  0.1    0.6  0.7
    - chr_col, pos_col: indicate column names of chromosome, position and SNP id, if not "CHR", "POS" and "SNP"
    - id_col: name of the ID column in the returned dosage dataframe. Default is 'genotype_id'
    
    Return:
    Return None if CHR and POS are not columns in the SNP file.
    Otherwise, return a Dataframe of dosage with sample IDs and SNP IDs labeled.
    Results have shape of sample x features(SNPs)
    
    # Path and file name examples
    vcf_path='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs'
    vcf_fn = 'max_unrelated_set_chr*.vcf.gz'
    snp_path = '/data100t1/home/wanying/CCHC/lipidomics/20231211_rerun/outputs/fastGWA/lipid_species_filter_by_pval_1e-05'
    snp_fn = 'TG-58:8-_[NL-22:6]_SNPs.pval_1e-05.txt'
    '''
    df_snps = pd.read_csv(snp_fn, sep='\t', dtype='str')
    
    # 'CHR' and 'POS' columns should be found in df_snps.columns
    if not chr_col in df_snps.columns:
        msg = '# Error: Chromosome column %s not found in the filtered SNP file' % chr_col
        logging.info(msg)
        return
    if not pos_col in df_snps.columns:
        msg = '# Error: Position column "%s" not found in the filtered SNP file' % pos_col
        logging.info(msg)
        return
    
    # Get list of positions from the SNP dataframe
    df_snps['position_for_tabix'] = 'chr' + df_snps[chr_col] + ':' + df_snps[pos_col] + '-' + df_snps[pos_col]
    dosage_all, snp_ids_all = [], [] # Store dosages and SNP ids of each chromosome
    df_snps[chr_col] = pd.to_numeric(df_snps[chr_col], downcast='integer') # Convert chromosome column to numeric value so that the process is ordered from chr1 to chr22 (or chr23)
    for i, (chr_num, df) in enumerate(df_snps.groupby(chr_col)):
        positions = df['position_for_tabix'].values
        if i==0:
            sample_ids, snp_ids, dosages = get_dosage_matrix(vcf_fn=os.path.join(vcf_path, vcf_fn.replace('*', str(chr_num))),
                                                             pos=positions, get_sample_ids=True, drop_multiallelic=drop_multiallelic)
        else:
            _, snp_ids, dosages = get_dosage_matrix(vcf_fn=os.path.join(vcf_path, vcf_fn.replace('*', str(chr_num))),
                                                             pos=positions, drop_multiallelic=drop_multiallelic)
        dosage_all += dosages
        snp_ids_all += snp_ids
    
    return pd.DataFrame(data=np.array(dosage_all),
                        columns=sample_ids,
                        index=snp_ids_all).T.reset_index().rename(columns={'index':'genotype_id'})