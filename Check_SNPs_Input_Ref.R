#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(TRUE)
inpfile <- args[1]
bimfile <- args[2]

inpdata <- data.table::fread(inpfile, header=T)
CN <- colnames(inpdata)
colnames(inpdata) <- c('rsid', 'chromosome', 'position', 'maf', 'beta', 'se')

bimdata <- data.table::fread(bimfile, header=F)
colnames(bimdata) <- c('chromosome', 'rsid', 'dummypos', 'position', 'allele1', 'allele2')

mergeDF <- dplyr::inner_join(inpdata, bimdata)
mergeDF <- mergeDF[, c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se')]
colnames(mergeDF) <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se')

## use space as the separator - recommended for FINEMAP
write.table(mergeDF, inpfile, row.names=F, col.names=T, sep=" ", quote=F, append=F)

