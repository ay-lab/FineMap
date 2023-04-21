#!/usr/bin/env Rscript

library(data.table)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(TRUE)
## input GWAS file
inpfile <- args[1]
## output file (converted GWAS)
outfile <- args[2]
## column containing the odd ratio information
ORcol <- as.integer(args[3])

## read the complete GWAS file
GWASData <- data.table::fread(inpfile, header=T)
cat(sprintf("\n Number of GWAS entries : %s ", nrow(GWASData)))

## estimate beta from odd ratio
## beta = ln(OR)
## check https://www.researchgate.net/post/How_can_I_calculate_beta_coefficient_and_its_error_from_Odds_Ratio_from_GWAS_summary_Statisitcs
GWASData$beta <- log(GWASData[, ORcol])
write.table(GWASData, outfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)



