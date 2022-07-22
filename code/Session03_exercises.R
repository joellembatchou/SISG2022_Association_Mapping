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

## Make a visual of the distribution of the phenotype.


# Question 3: Using PLINK, perform a GWAS using the phenotype file `sim_rels_pheno.txt` and the `sim_rels_geno.{bed,bim,fam}` genotype files.
system("plink2 --bfile sim_rels_geno --pheno sim_rels_pheno.txt --pheno-name <pheno_name> --maf <min_MAF> --geno <max_miss> --hwe <hwe_p_thresh> --glm allow-no-covars --out <output_prefix>")
  

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

# Question 6: Compute the genomic control inflation factor lambda_GC based on the p-values. Is there evidence of possible inflation due to confounding?


# Question 7: Now use REGENIE to perform a GWAS of the transferrin phenotype using a whole genome regression model.
## Write the list of samples with non-missing phenotype to file (FID/IID columns).

## Using PLINK, apply QC filters to remove variants with MAF below 5%, missingness above 1%, HWE p-value below 0.001, minor allele count (MAC) below 20.
### (hint: use `--write-snplist` to store list of variants passing QC without making a new BED file)

## Run REGENIE Step 1 to fit the null model and obtain polygenic predictions using a leave-one-chromosome-out (LOCO) scheme.
system("regenie --bed sim_rels_geno --phenoFile sim_rels_pheno.txt --step 1 --loocv --bsize 1000 --qt --extract <plink_QC_pass_snplist> --out <output_prefix_step1>")

## Run REGENIE Step 2 to perform association testing at the same variants analyzed in PLINK.
system("regenie --bed sim_rels_geno --phenoFile sim_rels_pheno.txt --step 2 --bsize 400 --qt  --pred <output_prefix_step1>_pred.list --extract <plink_GWAS_snplist> --out <output_prefix_step2>")

## Generate Manhatthan and Q-Q plots based on the association results and compute lambda_{GC}