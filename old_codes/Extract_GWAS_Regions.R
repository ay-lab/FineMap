#!/usr/bin/env Rscript

suppressMessages(library(GenomicRanges))
library(data.table)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {
	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(ov1[,1])
		ov_idx_file2 <- unique(ov1[,2])		
	} else {
		ov_idx_file1 <- ov1[,1]
		ov_idx_file2 <- ov1[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)
	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)
}

args <- commandArgs(TRUE)
## input GWAS file
inpfile <- args[1]
## output directory
outdir <- args[2]
## column containing the chromosome information
chrcol <- as.integer(args[3])
## column containing the SNP position information
poscol <- as.integer(args[4])
## column containing the p-value 
pvalcol <- as.integer(args[5])
## argument showing the span around a significant variant
## default = 3 Mb
if (length(args) > 5) {
	OFFSET <- as.integer(args[6])
} else {
	OFFSET <- 3000000
}

system(paste("mkdir -p", outdir))

regionfile <- paste0(outdir, '/GWAS_Regions.txt')
tempregionfile <- paste0(outdir, '/temp_GWAS_Regions.txt')

## read the complete GWAS file
GWASData <- data.table::fread(inpfile, header=T)
cat(sprintf("\n Number of GWAS entries : %s ", nrow(GWASData)))

## if the chromosome column do not contain the "chr" string (i.e. has chromosome numbers only), append
if (grepl("chr", GWASData[1, chrcol]) == FALSE) {
	GWASData[, chrcol] <- paste0("chr", GWASData[,chrcol])
}

## extract GWAS significant SNPs
GWASSigFile <- paste0(outdir, '/input_GWAS_SIG_pval_lt_5eMinus8.txt')
idx <- which(as.numeric(GWASData[, pvalcol]) < 0.00000005)
if (length(idx) > 0) {
	GWASSigData <- GWASData[idx, ]
	cat(sprintf("\n\n *** Number of GWAS significant SNPs (p-value < 5e-8) : %s \n\n", length(idx)))
	write.table(GWASSigData, GWASSigFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
} else {
	cat(sprintf("\n\n *** No GWAS significant SNPs (p-value < 5e-8) - exit !!! \n\n"))
	return()
}

## for each significant variant, extract 3 Mb region surrounding it
## and then merge those regions using bedtools
## each line of the generated file denotes individual non-overlapping regions
tempRegionData <- GWASSigData[, c(chrcol, poscol, poscol)]
colnames(tempRegionData) <- c('chr', 'start', 'end')
tempRegionData[, 2] <- tempRegionData[, 2] - (OFFSET / 2)
tempRegionData[, 3] <- tempRegionData[, 3] + (OFFSET / 2)
idx <- which(tempRegionData[, 2] < 0)
if (length(idx) > 0) {
	tempRegionData[idx, 2] <- 0
}
tempRegionData <- tempRegionData[order(tempRegionData[,1], tempRegionData[,2]), ]
write.table(tempRegionData, tempregionfile, row.names=F, col.names=F, sep="\t", quote=F, append=F)

## bedtools merge
system(paste0("bedtools merge -i ", tempregionfile, " > ", regionfile))

## now divide the non-overlapping regions and extract the summary statistics for individual regions
regiondata <- data.table::fread(regionfile, header=F)
GWAS_RegionDir <- paste0(outdir, '/GWAS_Regions')
system(paste("mkdir -p", GWAS_RegionDir))
for (i in 1:nrow(regiondata)) {
	currregiondata <- as.data.frame(regiondata[i, ])
	ov <- Overlap1D(GWASData[, c(chrcol, poscol, poscol)], currregiondata, boundary=0)
	curroutfile <- paste0(GWAS_RegionDir, '/Region_', i, '.txt')
	write.table(GWASData[ov$A_AND_B, ], curroutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
}

## remove temporary files
if (file.exists(tempregionfile)) {
	system(paste("rm", tempregionfile))
}



