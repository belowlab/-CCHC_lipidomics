'''
python xxx.py \
--train_vcf_path /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs \
--train_vcf_fn max_unrelated_set_chr*.vcf.gz \
--test_vcf_path /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_testing_vcfs \
--test_vcf_fn test_set_chr*.vcf.gz \
--snp_fn xxx \
--lipidomics_fn \
--lipid \
--output_path \
--output_prefix \
--model elastic_net
'''

import pandas as pd
import numpy as np
import os
import datetime
import subprocess
import argparse
import logging

from load_dosage import get_dosage_of_snps_from_gwas_output

# Train regression model, evaluate on test set

# ############### Process arguments ###############
parser = argparse.ArgumentParser(prog='Lipidomics modeling',
                                 description='Train and evaluate ML model with lipidomics data')

parser.add_argument('--train_vcf_path', type=str, default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs',
                    help='Path to the genotype vcf files for training')
parser.add_argument('--train_vcf_fn', type=str,
                    help='File name of the genotype vcf for training. Replace chromosome number with *. Eg. max_unrelated_set_chr*.vcf.gz')
parser.add_argument('--test_vcf_path', type=str, default='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_test_vcfs',
                    help='Path to the genotype vcf files for testing (Create predictions for the test data)')
parser.add_argument('--test_vcf_fn', type=str,
                    help='File name of the genotype vcf for testing. Replace chromosome number with *. Eg. test_set_chr*.vcf.gz')
parser.add_argument('--snp_fn', type=str,
                    help='A file contains SNPs to be loaded by the model. Must contain CHR and POS columns. Eg. TG-58:8-_[NL-22:6]_SNPs.pval_1e-05.txt') 

parser.add_argument('--lipidomics_fn', type=str,
                    help='A file contains lipidomics measures. May have different set of samples as long as there are some overlap with the training vcf.')
parser.add_argument('--lipid', type=str,
                    help='Name of the lipid to train and test model. Should be found in the file specified --lipidomics_fn')

parser.add_argument('--output_path', type=str, default='./', help='Path of the output files')
parser.add_argument('--output_prefix', type=str, help='Prefix of the output file')
parser.add_argument('--model', type=str, help='ML model type', choices=['elastic_net', 'ridge', 'adaB', 'gdB', 'mlp'])

args = parser.parse_args()


# #################### Get dosages ####################
df_dosage = get_dosage_of_snps_from_gwas_output()

# ########## Load lipidomic measures and match order with the genotype dosages ##########
df_lipid = pd.read_csv(aargs.lipidomics_fn, sep='\t')
# Assume the first column is the ID column
df_lipid = df_lipid.set_index(keys=df_lipid.columns[0]).reindex(index=df_dosage.index)
