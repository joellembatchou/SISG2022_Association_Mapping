require(data.table)
require(dplyr)
require(tidyr)
require(GWASTools)
require(ggplot2)

####################################################
## GWAS in Samples with Structure & Using REGENIE ##
####################################################
# Question 1. Examine the dataset
## How many samples are present? 

## How many SNPs? In how many chromosomes?


# Question 2. Examine the phenotype data:
## How many individuals in the study have measurements?

## What is the distribution of the phenotype? (hint: plot a histogram)


# Question 3: Using PLINK, perform a GWAS of transferrin using the phenotype file `Transferrin_pheno.txt` and the `Transferrin_height.{bed,bim,fam}` genotype files. Only perform association test on SNPs that pass the following quality control threshold filters:
system("plink2 --bfile Transferrin_height --pheno Transferrin_pheno.txt --pheno-name <pheno_name> --maf <min_MAF> --geno <max_miss> --hwe <hwe_p_thresh> --glm allow-no-covars --out <output_prefix>")
  

# Question 4: Make a Manhattan plot of the association results using the `manhattanPlot()` R function
manhattanPlot(
  p = <pvalues>,
  chromosome = <chromosomes>, 
  thinThreshold = 1e-4,
  main= <title>
)

# Question 5: Make a Q-Q plot of the association results using the `qqPlot()` R function
qqPlot(
  pval = <pvalues>,
  thinThreshold = 1e-4,
  main= <title>
)

# Question 6: Compute the genomic control inflation factor $\lambda_{GC}$ based on the p-values.


# Question 7: Now use REGENIE to perform a GWAS of the transferrin phenotype using a whole genome regression model.
## Write the list of samples with non-missing phenotype to file (FID/IID columns).

## Using PLINK, apply QC filters to remove variants with MAF below 5%, missingness above 1%, HWE p-value below 0.001, minor allele count (MAC) below 20. Make sure to specify the ID of samples to analyze using `--keep`. (hint: use `--write-snplist` to store list of variants passing QC without making a new BED file)

## Download the Step 1 output file

## Generate the LOCO list file

## Run REGENIE Step 2 to perform association testing at the same variants analyzed in PLINK

## Generate Manhatthan and Q-Q plots based on the association results and compute lambda_{GC}