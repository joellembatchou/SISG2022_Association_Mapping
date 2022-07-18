require(data.table)
require(dplyr)
require(tidyr)
require(bigsnpr)
require(ggplot2)

####################################
## Population Structure Inference ##
####################################
# Question 1. Examine the dataset
## How many samples are present? 

## How many SNPs?

## What is the number of samples in each population?


# Question 2. Get the first 10 principal components (PCs) in PLINK using all SNPs.
system("plink --bfile <plink_bed_prefix> --pca 10 --out <output_prefix>")

## Make a scatterplot of the first two PCs with each point colored by population membership. 

## Interpret the first two PCs, what ancestries are they reflecting?

## Make a scree plot of the eigenvalues for the first 10 PCs. Approximate the proportion of variance explained by the first two PCs.


# Question 3: Now redo Question 2 above using the bigsnpr R package specifying a R2 threshold of 0.2 (i.e. LD pruning) as well as a minimum minor allele count (MAC) of 20.
# Hint: Use the relevel function to create a new genotype vector with reference genotype TT.
obj.bed <- bed(bedfile = <plink_bed_file>)
pc.out <- bed_autoSVD(
  obj.bed, 
  thr.r2 = <r2_threshold>, 
  k = <number_of_PCs>, 
  min.mac = <min_MAC>
)

## Make a scatter plot of the first two principal components (PCs) with each point colored according to population membership. 
plot(pc.out, type = "scores", scores = 1:2) + aes(color = ...)

## Does the plot change from the one in Question 2?

## Check the SNP loadings for the first 10 PCs.
plot(pc.out, type = "loadings", scores = 1:<number_of_PCs>, coeff = 0.4)


# Question 4: Predict the proportional Native American and European Ancestry for the HapMap MXL from the PCA output in Question 3 *using one of the principal components*. (Which PC is most appropriate for this analysis?) Assume that the HapMap MXL have negligible African Ancestry.


# Question 5: Make a barplot of the proportional ancestry estimates from question 4.
# hint: your data frame should have sample index, population (CEU/NAM) and the corresponding ancestry proportion estimate
# hint: make a stacked barplot using `geom_bar(position="stack", stat="identity")`


# Extra - Question 6: Check if there are samples related 2nd degree or closer. If so, run PCA as in Question 3 removing these samples then project the remaining samples onto the PC space.
# identify 2nd degree relateds or closer
rel.df <- snp_plinkKINGQC(
  plink2.path = <plink2_path>, 
  bedfile.in = <plink_bed_prefix>, 
  thr.king = 2^-3.5,
  make.bed = FALSE
)
# run pca only on unrelateds
pc.out.unrels <- bed_autoSVD(
  obj.bed, 
  ind.row = <...>,
  thr.r2 = <r2_threshold>, 
  k = <number_of_PCs>, 
  min.mac = <min_MAC>
)
# project relateds onto PC space
proj.rels <- bed_projectSelfPCA(
  pc.out.unrels, 
  obj.bed,
  ind.row = <...>
)
# `proj.rels$OADP_proj` will contain PC projections for related samples