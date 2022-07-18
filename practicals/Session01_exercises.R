require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)

######################################
## Case-Control Association Testing ##
######################################
# Read in the data
LHON.df <- fread("https://raw.githubusercontent.com/joellembatchou/SISG2022_Association_Mapping/master/data/LHON.txt", header=TRUE)

# Question 1. Examine the variables in the dataset:
## How many observations? 

## How many cases/controls?

## What is the distribution of the genotypes across cases/controls?

## What about for allele types?


# Question 2. Perform a logistic regression analysis for this data with `CC` as the reference genotype using the `glm()` function.


# Question 3: Redo the logistic regression analysis from question 2, but with `TT` as the reference genotype.
# Hint: Use the relevel function to create a new genotype vector with reference genotype TT.


# Question 4: Is there evidence of differences in odds of being a case for the `CT` and `TT` genotypes (compared to `CC`)?


##################################################
## Association Testing with Quantitative Traits ##
##################################################
# Read in the data
BP.df <- fread("https://raw.githubusercontent.com/joellembatchou/SISG2022_Association_Mapping/master/data/bpdata.csv", header=TRUE)

# Question 1. Perform a linear regression of systolic blood pressure (`sbp`) on `SNP3` using the `lm()` function.
## Additive (linear) model 

## Dominant model 

## Recessive model

## 2 parameter model


# Question 2. Provide a plot illustrating the relationship between sbp and the three genotypes at SNP3.

# Question 3: Redo the linear regression analysis of `sbp` but this time adjust for `sex` and `bmi`.

# Question 4: What proportion of the heritability of `sbp` is explained by all of the 11 SNPs together?