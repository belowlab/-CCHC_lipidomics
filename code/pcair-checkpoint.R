library(SNPRelate)
library(GENESIS)
library(GWASTools)
# Refer to code here: https://www.bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html#plink-files
# HH's code is the same: /vgipiper04/CCHC/GENESIS/CCHC_GENESIS.log
# pcair documents: https://bioconductor.org/packages/3.17/bioc/manuals/GENESIS/man/GENESIS.pdf

# Read in genotype data
snpgdsBED2GDS(bed.fn = "/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_preimpute/lipidomic_subset.bed", 
              bim.fn = "/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_preimpute/lipidomic_subset.bim", 
              fam.fn = "/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_plink_preimpute/lipidomic_subset.fam", 
              out.gdsfn = "/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_pca/genotype.gds")

# ################ 1. LD pruning ################
# read in GDS data
gdsfile="/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_pca/genotype.gds"
gds <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1), verbose=TRUE, num.thread=64) # Seems multithreding is not supported, can ignore num.thread
pruned <- unlist(snpset, use.names=FALSE)
print(length(pruned))

# ################ 2. Pairwise Measures of Ancestry Divergence ################
# Run KING-robust: relationship inference in the presence of population stratification robust relationship inference across family
ibd.robust <- snpgdsIBDKING(gds, verbose=TRUE)
print(names(ibd.robust))
KINGmat <- kingToMatrix(ibd.robust)
KINGmat[1:5,1:5]

# ################ 3. Running PC-AiR ################
# run PC-AiR on pruned SNPs
snpgdsClose(gds) # Close the file so GdsGenotypeReader can open it
CCHC_geno <- GdsGenotypeReader(filename = gdsfile)
CCHC_genoData <- GenotypeData(CCHC_geno)
CCHC_genoData
mypcair <- pcair(CCHC_genoData, kinobj = KINGmat, divobj = KINGmat,
                 snp.include = pruned)
summary(mypcair)

# ################ 4. Output from PC-AiR ################
output_dir = "/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_pca/"
# vectors: A matrix of principal components; each column is a principal component
write.table(mypcair$vectors, paste(output_dir, 'CCHC_PCair.vectors', sep=''), col.names=F, row.names=T, sep='\t', quote=F)

# values: A vector of eigenvalues matching the principal components.
# These values are determined from the standard PCA run on the ’unrelated subset’
write.table(mypcair$values, paste(output_dir, 'CCHC_PCair.values.txt', sep=''), col.names=F, row.names=T, sep='\t', quote=F)

# unrels: A vector of IDs for individuals in the ’unrelated subset’
write.table(mypcair$unrels, paste(output_dir, 'CCHC_PCair.unel.txt', sep=''), col.names=F, row.names=F, sep='\t', quote=F)