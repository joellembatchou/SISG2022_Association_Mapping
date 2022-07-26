require(data.table)
require(dplyr)
require(tidyr)
require(GWASTools)
require(ggplot2)
require(patchwork)

## Set to your home directory if on server
setwd("~/")

## Data sets are in : /data/SISG2022M15/data/

####################################
## Case-control imbalance in GWAS ##
####################################

#### Data simulation - genotypes
N <- 10e3
# Generate a configuration file specifying allele frequencies (a,b) for Uniform(a,b) distribution
write(paste0("5000 common 0.05 0.5 1 1"), "sim.config")
write(paste0("5000 rare 0.001 0.01 1 1"), "sim.config", append = TRUE)
# Run PLINK1.9
system(paste0("/data/SISG2022M15/exe/plink1.9 --make-bed --simulate sim.config --simulate-ncases ", N, " --simulate-ncontrols 0 --simulate-prevalence 0.1  --out cc_imb_geno"))

#### Data simulation - phenotypes
# get FID/IID from FAM file
sample.ids <- fread("cc_imb_geno.fam", header = FALSE)
## Set prevalence = 10% (CCR 1:9)
y1 <- rbinom(N, 1, prob = 0.1 )
## Set prevalence = 1% (CCR 1:99)
y2 <- rbinom(N, 1, prob = 0.01 )
## Set prevalence = 0.5% (CCR 1:199)
y3 <- rbinom(N, 1, prob = 0.005 )
# write to file
data.frame(FID = sample.ids$V1, IID = sample.ids$V2, Y1 = y1, Y2 = y2, Y3 = y3) %>%
  fwrite("cc_imb_pheno.txt", sep = "\t", na = NA, quote = FALSE)



# Question 1. Run GWAS in REGENIE (step 2 only) analyzing all 3 traits as quantitative.
system("/data/SISG2022M15/exe/regenie --bed cc_imb_geno --phenoFile cc_imb_pheno.txt --step 2 --bsize 400 --qt --ignore-pred --out <output_prefix>")


# Question 2. Read in the three summary statistics files in R and make a QQ plot of the p-values for each phenotype. Since these are null SNPs, how does it compare to what we expect?
sumstats.y1 <- fread("<output_prefix>_Y1.regenie") %>% 
  mutate(Pval = 10^(-LOG10P), Z = sign(BETA) * sqrt(CHISQ)) # get Pvalues & Z-scores
qqPlot( pval = sumstats.y1$Pval )
# do for remaining traits as well (Y2, Y3)


# Question 3: Make a histogram of the test statistics for each phenotype and overlay with a normal distribution. How well do they match? What do you observe as the case-control imbalance gets more severe?

## R function to easily make this plot for different phenotypes.
plot.sumstats.hist <- function(df, title = ""){
  df %>%
    ggplot( aes(x = Z) ) +
    geom_histogram(aes(y = ..density..), colour="black", fill="white", bins = 100) +
    stat_function(
      fun = dnorm, 
      col = "red",
      args = list(mean = mean(df$Z), sd = sd(df$Z))
    ) +
    labs(title = title)
}
plot.sumstats.hist(sumstats.y1, title = "Y1")
# do for remaining traits as well (Y2, Y3)


# Question 4: Re-do 3 but now separate the histogram for common and rare SNPs. 
## First separate the data frame based on common/rare simulated SNPs.
sumstats.y1.common <- sumstats.y1[ grepl("common", ID), ]
sumstats.y1.rare <- sumstats.y1[ grepl("rare", ID), ]
# do for remaining traits as well (Y2, Y3)

## Make a histogram of the test statistics distribution at common/rare SNPs. What do you observe across the different case-control imbalances? 
plot.sumstats.hist(sumstats.y1.common, title = "Y1 - Common SNPs") | plot.sumstats.hist(sumstats.y1.rare, title = "Y1 - Rare SNPs")
