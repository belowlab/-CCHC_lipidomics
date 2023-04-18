# Generate IBD files in plink to get max unrelated set from primus
bfile=/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_preimpute/lipidomic_subset
plink --bfile ${bfile} --genome

# run primus
/local-data100t1/gapps/PRIMUS_v1.9.0/bin/run_PRIMUS.pl -p plink.genome 
