require(data.table)
require(dplyr)
require(tidyr)
require(BEDMatrix)
require(SKAT)
require(ACAT)
require(ggplot2)

## Set to your home directory if on server
setwd("~/")

## Data sets are in : /data/SISG2022M15/data/

###########################
## Rare Variant Analysis ##
###########################
# Question 1. Using PLINK, extract **rare variants** in a new PLINK BED file.
# PLINK BED file is in /data/SISG2022M15/data/rv_geno_chr1.bed
system("/data/SISG2022M15/exe/plink2 --bfile /data/SISG2022M15/data/rv_geno_chr1 --max-maf <..> --maj-ref force --make-bed --out <output_prefix>")


# Question 2. Load the data in R
## Read in the SNPs from the BED file created in Q1 using R function `BEDMatrix()
G <- BEDMatrix( "<bed_Q1_prefix>", simple_names = TRUE) %>% as.matrix
iid.geno <- rownames(G) # sample IDs of individuals in the genotype data

## Load the phenotype data from `rv_pheno.txt`
y <- fread("<pheno_file>", header = TRUE)

## Keep only samples who are present both in the genotype as well as phenotype data and who don't have missing values for the phenotype

# Question 3: Examine the genotype data:
## Compute the minor allele frequency (MAF) for each SNP and plot histogram. (hint: use `na.rm=TRUE` when calling `mean()`)
maf <- colMeans(G, na.rm = TRUE)/2

## Check for missing values


# Question 4: Run the single variant association tests in PLINK (only for the extracted variants).
system("/data/SISG2022M15/exe/plink2 --bfile <BED_file_with_extracted_SNPs> --pheno /data/SISG2022M15/data/rv_pheno.txt --pheno-name <pheno_name> --glm allow-no-covars --out <output_prefix>")
## What would be your significance threshold after applying Bonferroni correction for the multiple tests (assume the significance level is 0.05)? Is anything significant after this correction?

## Make a volcano plot (i.e. log10 p-values (y) vs effect sizes (x) ).


# Question 5: We will first compare three collapsing/burden approaches:
## CAST (Binary collapsing approach): for each individual, count where they have a rare allele at any of the sites
burden.cast <- ..
lm( <y> ~ burden.cast) %>% summary

## MZ Test/GRANVIL (Count based collapsing): for each individual, count the total number of sites where a rare allele is present
burden.mz <- ..
lm( <y> ~ burden.mz) %>% summary

## Weighted burden test: for each individual, take a weighted count of the rare alleles across sites
weights <- dbeta( maf , 1, 25)
burden.weighted <- ..
lm( <y> ~ burden.weighted) %>% summary


# Question 6: Now use SKAT to test for an association. 
skat.null <- SKAT_Null_Model( <phenotype_vector> ~ 1 , out_type = "C")
SKAT( <genotype_matrix>, skat.null )


# Question 7: Run the omnibus SKAT, but consider setting $\rho$ (i.e.`r.corr`) to 0 and then 1.
SKAT( <genotype_matrix>, skat.null, r.corr = <rho_value>)


# Question 8: Now the omnibus version of SKAT, but use the “optimal.adj” approach which searches across a range of rho values.
SKAT( <genotype_matrix>, skat.null, method="optimal.adj")


# Question 9: Run ACATV on the single variant p-values.
acat.weights <- weights * weights * maf * (1 - maf)
ACAT( <pvalues>, weights = acat.weights)


# Question 10: Run ACATO combining the SKAT and BURDEN p-values (from Question 7) with the ACATV p-value (from Question 9).
ACAT( c(<pvalue_SKAT>, <pvalue_Burden>, <pvalue_ACATV>) )

