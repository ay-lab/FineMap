#!/usr/bin/env Rscript

#========================
## this script is for fine-mapping eQTL / GWAS summary statistics
## a wrapper of the FINEMAP package
#========================
suppressMessages(library(data.table))
suppressMessages(library(yaml))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))

options(scipen = 10)
options(datatable.fread.datatable=FALSE) 

## significance threshold for GWAS
GWAS_SIG_THR <- 5e-8

##===========
## dummy SE value (small floating point value) for zero entries
DUMMY_SE_VAL <- 1e-10

## dummy p-value (for zero entries)
DUMMY_P_VAL <- 1e-50

##=================
## overlapping 1D genomic regions
##=================
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

##=============
## function to extract GWAS loci
##=============
Extract_GWAS_Regions <- function(GWASData, OFFSET, regionfile) {
	
	tempregionfile <- paste0(dirname(regionfile), '/temp_GWAS_Regions.txt')

	## extract GWAS significant SNPs	
	idx <- which(as.numeric(GWASData$pval_GWAS) < GWAS_SIG_THR)
	if (length(idx) > 0) {
		GWASSigData <- GWASData[idx, ]
	} else {
		cat(sprintf("\n\n *** No GWAS significant SNPs (p-value < 5e-8) - exit !!! \n\n"))
		return()
	}

	## for each significant variant, extract 3 Mb region surrounding it
	## and then merge those regions using bedtools
	## each line of the generated file denotes individual non-overlapping regions
	tempRegionData <- GWASSigData[, c(1, 2, 2)]
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
	GWAS_RegionDir <- paste0(dirname(regionfile), '/GWAS_Regions')
	system(paste("mkdir -p", GWAS_RegionDir))
	for (i in 1:nrow(regiondata)) {
		currregiondata <- as.data.frame(regiondata[i, ])
		ov <- Overlap1D(GWASData[, c(1, 2, 2)], currregiondata, boundary=0)
		curroutfile <- paste0(GWAS_RegionDir, '/Region_', i, '.txt')
		write.table(GWASData[ov$A_AND_B, ], curroutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}

	## remove temporary files
	if (file.exists(tempregionfile)) {
		system(paste("rm", tempregionfile))
	}
}

##==============
## function to filter SNPs according to the input genotype file
## used together with LDStore
##==============
Filter_LDStore_SNP <- function(inpfile, bimfile) {
	inpdata <- data.table::fread(inpfile, header=T)
	bimdata <- data.table::fread(bimfile, header=F)
	colnames(bimdata) <- c('chromosome', 'rsid', 'dummypos', 'position', 'allele1', 'allele2')
	mergeDF <- dplyr::inner_join(inpdata, bimdata)
	mergeDF <- mergeDF[, c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se')]
	colnames(mergeDF) <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se')
	## use space as the separator - recommended for FINEMAP
	write.table(mergeDF, inpfile, row.names=F, col.names=T, sep=" ", quote=F, append=F)
}

##==============
## function to compute LD from GWAS SNPs
## with respect to individual loci
## and then run FINEMAP
##==============
Run_Finemap_GWAS <- function(FINEMAPInpDir, BaseOutDir_GWAS, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, GENOTYPEDIR, ldstoreexec, finemapexec, samplecount, NUMCAUSALSNP, NUMTHREAD) {	
		
	##=== process individual GWAS regions
	i <- 0 	
	while(1) {
		i <- i + 1
		gwasfile <- paste0(BaseOutDir_GWAS, '/GWAS_Regions/Region_', i, '.txt')
		if (file.exists(gwasfile) == FALSE) {
			break
		}
		currgwasdata <- data.table::fread(gwasfile, header=T)
		cat(sprintf("\n\n==>> processing gwasfile : %s Number of GWAS entries : %s ", gwasfile, nrow(currgwasdata)))

		finemap_z_file <- paste0(FINEMAPInpDir, '/Region', i, '.z')

		## master file which will be used as the input of FINEMAP
		Masterfile <- paste0(FINEMAPInpDir, '/Region', i, '.master')
		currstr <- "z;bgen;bgi;bcor;ld;n_samples;bdose"
		write(currstr, file=Masterfile, append=FALSE)

		## get the chromosome
		currchr <- currgwasdata[1, 1]
		## get the chromosome number
		currchrnum <- currchr
		currchrnum <- as.integer(gsub("chr", "", currchrnum))
		cat(sprintf("\n current chromosome : %s number : %s ", currchr, currchrnum))
		
		## reference 1000G genotypes - variants are referred by chrNum:pos
		## we are using all the reference panels
		bgenfile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bgen')
		bgenbgifile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bgen.bgi')
		bimfile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bim')
		if ((file.exists(bgenfile) == FALSE) | (file.exists(bgenbgifile) == FALSE) | (file.exists(bimfile) == FALSE)) {
			next
		}		

		## use chromosome:position format for variant ID
		## this file will also be useful for FINEMAP
		## Note: output file is a space separated file
		## Note: if allele frequency column (afcol) is provided, insert MAF (< 0.5) and convert the allele frequency if required
		## otherwise, use a dummy value like 0.1
		## also check the standard error - if 0, insert a small floating point positive value
		if ("rsID" %in% colnames(currgwasdata)) {
			finemapdf <- data.frame(rsid=currgwasdata$rsID, chromosome=rep(currchrnum, nrow(currgwasdata)), position=currgwasdata$pos)
		} else {
			finemapdf <- data.frame(chromosome=rep(currchrnum, nrow(currgwasdata)), position=currgwasdata$pos)
			finemapdf$rsid <- paste(finemapdf$chromosome, finemapdf$position, sep=":")
			finemapdf <- finemapdf[, c("rsid", "chromosome", "position")]
		}
		finemapdf$allele1 <- currgwasdata$allele1
		finemapdf$allele2 <- currgwasdata$allele2
		if ("AF_GWAS" %in% colnames(currgwasdata)) {
			finemapdf$maf <- pmin(currgwasdata$AF_GWAS, (1 - currgwasdata$AF_GWAS))
		} else {
			finemapdf$maf <- rep(0.1, nrow(currgwasdata))
		}
		finemapdf$beta <- currgwasdata$beta_GWAS
		finemapdf$se <- currgwasdata$SE_gwas
		write.table(finemapdf, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F)

		##========= before running LDStore, check if these SNPs are present in the reference genotype (PLINK) output
		##========= otherwise filter those entries
		n <- as.integer(system(paste("cat", finemap_z_file, "| wc -l"), intern = TRUE))		
		cat(sprintf("\n before filtering SNPs from FINEMAP input - number of entries : %s ", n))

		Filter_LDStore_SNP(finemap_z_file, bimfile)
		n <- as.integer(system(paste("cat", finemap_z_file, "| wc -l"), intern = TRUE))
		cat(sprintf("\n after filtering SNPs from FINEMAP input - number of entries : %s ", n))

		##======== another filtering - sourya
		##======== sometimes same SNPs (same rsID, position) exist multiple times, even in the original GWAS data - filter them
		finemapdf <- data.table::fread(finemap_z_file, header=T)
		finemapdf2 <- finemapdf[!duplicated(finemapdf$rsid), ]
		cat(sprintf("\n after filtering SNPs from FINEMAP input - number of entries : %s ", nrow(finemapdf2)))
		write.table(finemapdf2, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F)

		##===============
		## run LDStore
		##===============

		##========= write the master file entries
		outbcorfile <- paste0(FINEMAPInpDir, '/Region', i, '.bcor')
		finemap_ld_file <- paste0(FINEMAPInpDir, '/Region', i, '.ld')
		bdosefile <- paste0(FINEMAPInpDir, '/Region', i, '.bdose')
		currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount, ";", bdosefile)
		write(currstr, file=Masterfile, append=TRUE)

		## now execute ldstore using the generared master file
		system(paste(ldstoreexec, "--in-files", Masterfile, "--read-only-bgen --write-text --write-bdose --n-threads", NUMTHREAD, "--memory 40"))

		##===============
		## run Fine mapping
		##===============

		## master file containing the finemap execution input
		## when --cond (stepwise conditioning) option is used 
		Masterfile_cond <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.master')
		currstr <- "z;ld;snp;config;cred;log;n_samples"
		write(currstr, file=Masterfile_cond, append=FALSE)

		## master file containing the finemap execution input
		## when --sss (shotgun stochastic search) option is used 
		Masterfile_sss <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.master')
		currstr <- "z;ld;snp;config;cred;log;n_samples"
		write(currstr, file=Masterfile_sss, append=FALSE)

		## construct the master file "Masterfile_cond"
		finemap_out_snp_file <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.snp')
		finemap_out_config_file <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.config')
		finemap_out_cred_file <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.cred')
		finemap_out_log_file <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.log')
		currstr <- paste0(finemap_z_file, ";", finemap_ld_file, ";", finemap_out_snp_file, ";", finemap_out_config_file, ";", finemap_out_cred_file, ";", finemap_out_log_file, ";", samplecount)
		write(currstr, file=Masterfile_cond, append=TRUE)
		
		## construct the master file "Masterfile_sss"
		finemap_out_snp_file <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.snp')
		finemap_out_config_file <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.config')
		finemap_out_cred_file <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.cred')
		finemap_out_log_file <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.log')
		currstr <- paste0(finemap_z_file, ";", finemap_ld_file, ";", finemap_out_snp_file, ";", finemap_out_config_file, ";", finemap_out_cred_file, ";", finemap_out_log_file, ";", samplecount)
		write(currstr, file=Masterfile_sss, append=TRUE)

		## fine mapping when --cond (stepwise conditioning) option is used 
		system(paste(finemapexec, "--cond --in-files", Masterfile_cond, "--log --n-causal-snps", NUMCAUSALSNP, "--n-threads", NUMTHREAD))

		## fine mapping when --sss (shotgun stochastic search) option is used 
		system(paste(finemapexec, "--sss --in-files", Masterfile_sss, "--log --n-causal-snps", NUMCAUSALSNP, "--n-threads", NUMTHREAD))

	}	# end GWAS loci loop
}

##==============
## function to run Fine mapping on eQTL data, for individual genes
##==============
Run_Finemap_eQTL <- function(eQTLData, OutDir, FINEMAPInpDir, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, samplecount, NUMTHREAD, GENOTYPEDIR, ldstoreexec, finemapexec) {

	## process for individual genes
	if ("GeneName" %in% colnames(eQTLData)) {
		GeneList <- unique(eQTLData$GeneName)
		GeneCol <- which(colnames(eQTLData) == "GeneName")
	} else {
		GeneList <- unique(eQTLData$GeneID)
		GeneCol <- which(colnames(eQTLData) == "GeneID")
	}
	cat(sprintf("\n\n ==>> Number of entries in eQTL data : %s Number of unique genes in eQTL data : %s  column containing gene info : %s ", nrow(eQTLData), length(GeneList), GeneCol))

	for (i in 1:length(GeneList)) {
		currgene <- GeneList[i]
		currgene_eQTLdata <- eQTLData[which(eQTLData[, GeneCol] == currgene), ]
		cat(sprintf("\n eQTL data - processing gene : %s  number of entries in eQTL data for this gene : %s ", currgene, nrow(currgene_eQTLdata)))

		FINEMAPInpDir_currgene <- paste0(FINEMAPInpDir, '/', currgene)
		system(paste("mkdir -p", FINEMAPInpDir_currgene))

		## master file which will be used as the input of FINEMAP
		finemap_z_file <- paste0(FINEMAPInpDir_currgene, '/Region', i, '.z')
		Masterfile <- paste0(FINEMAPInpDir_currgene, '/Region', i, '.master')
		currstr <- "z;bgen;bgi;bcor;ld;n_samples;bdose"
		write(currstr, file=Masterfile, append=FALSE)

		## get the chromosome
		currchr <- currgene_eQTLdata[1, 1]
		## get the chromosome number
		currchrnum <- currchr
		currchrnum <- as.integer(gsub("chr", "", currchrnum))
		cat(sprintf("\n current chromosome : %s number : %s ", currchr, currchrnum))

		## reference 1000G genotypes - variants are referred by chrNum:pos
		## we are using all the reference panels
		bgenfile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bgen')
		bgenbgifile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bgen.bgi')
		bimfile <- paste0(GENOTYPEDIR, '/all_phase3_1000G_', currchr, '.bim')
		if ((file.exists(bgenfile) == FALSE) | (file.exists(bgenbgifile) == FALSE) | (file.exists(bimfile) == FALSE)) {
			next
		}		

		## use chromosome:position format for variant ID
		## this file will also be useful for FINEMAP
		## Note: output file is a space separated file
		## Note: if allele frequency column (afcol) is provided, insert MAF (< 0.5) and convert the allele frequency if required
		## otherwise, use a dummy value like 0.1
		## also check the standard error - if 0, insert a small floating point positive value
		if ("rsID" %in% colnames(currgene_eQTLdata)) {
			finemapdf <- data.frame(snpid=currgene_eQTLdata$rsID, chromosome=rep(currchrnum, nrow(currgene_eQTLdata)), position=currgene_eQTLdata$pos)
		} else {
			finemapdf <- data.frame(chromosome=rep(currchrnum, nrow(currgene_eQTLdata)), position=currgene_eQTLdata$pos)
		}
		finemapdf$rsid <- paste(finemapdf$chromosome, finemapdf$position, sep=":")
		finemapdf$chromosome <- as.integer(finemapdf$chromosome)
		finemapdf <- finemapdf[, c("rsid", "chromosome", "position")]		
		if ("AF_eQTL" %in% colnames(currgene_eQTLdata)) {
			finemapdf$maf <- pmin(currgene_eQTLdata$AF_eQTL, (1 - currgene_eQTLdata$AF_eQTL))
		} else {
			finemapdf$maf <- rep(0.1, nrow(currgene_eQTLdata))
		}
		finemapdf$beta <- currgene_eQTLdata$beta_eQTL
		finemapdf$se <- currgene_eQTLdata$SE_eQTL
		write.table(finemapdf, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F)

		##========= before running LDStore, check if these SNPs are present in the reference genotype (PLINK) output
		##========= otherwise filter those entries
		n <- as.integer(system(paste("cat", finemap_z_file, "| wc -l"), intern = TRUE))		
		cat(sprintf("\n before filtering SNPs from FINEMAP input - number of entries : %s ", n))

		Filter_LDStore_SNP(finemap_z_file, bimfile)
		n <- as.integer(system(paste("cat", finemap_z_file, "| wc -l"), intern = TRUE))
		cat(sprintf("\n after filtering SNPs from FINEMAP input - number of entries : %s ", n))

		##======== another filtering - sourya
		##======== sometimes same SNPs (same rsID, position) exist multiple times, even in the original GWAS data - filter them
		finemapdf <- data.table::fread(finemap_z_file, header=T)
		finemapdf2 <- finemapdf[!duplicated(finemapdf$rsid), ]
		cat(sprintf("\n after filtering SNPs from FINEMAP input - number of entries : %s ", nrow(finemapdf2)))
		write.table(finemapdf2, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F)

		##===============
		## run LDStore
		##===============

		##========= write the master file entries
		outbcorfile <- paste0(FINEMAPInpDir_currgene, '/Region', i, '.bcor')
		finemap_ld_file <- paste0(FINEMAPInpDir_currgene, '/Region', i, '.ld')
		bdosefile <- paste0(FINEMAPInpDir_currgene, '/Region', i, '.bdose')
		currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount, ";", bdosefile)
		write(currstr, file=Masterfile, append=TRUE)

		## now execute ldstore using the generared master file
		system(paste(ldstoreexec, "--in-files", Masterfile, "--read-only-bgen --write-text --write-bdose --n-threads", NUMTHREAD, "--memory 40"))

		##===============
		## run Finemap
		##===============

		FINEMAPOutDir_cond_currgene <- paste0(BaseOutDir_FINEMAP_cond, '/', currgene)
		system(paste("mkdir -p", FINEMAPOutDir_cond_currgene))

		FINEMAPOutDir_sss_currgene <- paste0(BaseOutDir_FINEMAP_sss, '/', currgene)
		system(paste("mkdir -p", FINEMAPOutDir_sss_currgene))

		## master file containing the finemap execution input
		## when --cond (stepwise conditioning) option is used 
		Masterfile_cond <- paste0(FINEMAPOutDir_cond_currgene, '/Region', i, '.master')
		currstr <- "z;ld;snp;config;cred;log;n_samples"
		write(currstr, file=Masterfile_cond, append=FALSE)

		## master file containing the finemap execution input
		## when --sss (shotgun stochastic search) option is used 
		Masterfile_sss <- paste0(FINEMAPOutDir_sss_currgene, '/Region', i, '.master')
		currstr <- "z;ld;snp;config;cred;log;n_samples"
		write(currstr, file=Masterfile_sss, append=FALSE)

		## construct the master file "Masterfile_cond"
		finemap_out_snp_file <- paste0(FINEMAPOutDir_cond_currgene, '/Region', i, '.snp')
		finemap_out_config_file <- paste0(FINEMAPOutDir_cond_currgene, '/Region', i, '.config')
		finemap_out_cred_file <- paste0(FINEMAPOutDir_cond_currgene, '/Region', i, '.cred')
		finemap_out_log_file <- paste0(FINEMAPOutDir_cond_currgene, '/Region', i, '.log')
		currstr <- paste0(finemap_z_file, ";", finemap_ld_file, ";", finemap_out_snp_file, ";", finemap_out_config_file, ";", finemap_out_cred_file, ";", finemap_out_log_file, ";", samplecount)
		write(currstr, file=Masterfile_cond, append=TRUE)

		## construct the master file "Masterfile_sss"
		finemap_out_snp_file <- paste0(FINEMAPOutDir_sss_currgene, '/Region', i, '.snp')
		finemap_out_config_file <- paste0(FINEMAPOutDir_sss_currgene, '/Region', i, '.config')
		finemap_out_cred_file <- paste0(FINEMAPOutDir_sss_currgene, '/Region', i, '.cred')
		finemap_out_log_file <- paste0(FINEMAPOutDir_sss_currgene, '/Region', i, '.log')
		currstr <- paste0(finemap_z_file, ";", finemap_ld_file, ";", finemap_out_snp_file, ";", finemap_out_config_file, ";", finemap_out_cred_file, ";", finemap_out_log_file, ";", samplecount)
		write(currstr, file=Masterfile_sss, append=TRUE)

		## fine mapping when --cond (stepwise conditioning) option is used 
		system(paste(finemapexec, "--cond --in-files", Masterfile_cond, "--log --n-causal-snps", NUMCAUSALSNP, "--n-threads", NUMTHREAD))

		## fine mapping when --sss (shotgun stochastic search) option is used 
		system(paste(finemapexec, "--sss --in-files", Masterfile_sss, "--log --n-causal-snps", NUMCAUSALSNP, "--n-threads", NUMTHREAD))

	}	# end gene loop
}

##===================
## function - part of summarization - writes fine-mapped final summary SNPs
##===================

Consolidate_Summary_Output <- function(outdir) {
	system(paste0("cat ", outdir, "/Region_*/snp_cum_prob_95.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", outdir, "/FINAL_snp_cum_prob_95.txt"))
	system(paste0("cat ", outdir, "/Region_*/summary_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", outdir, "/FINAL_summary_credible_set.txt"))
	system(paste0("cat ", outdir, "/Region_*/top_snp_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", outdir, "/FINAL_top_snp_credible_set.txt"))
	system(paste0("cat ", outdir, "/Region_*/top_snp_log10BF_geq_2.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", outdir, "/FINAL_top_snp_log10BF_geq_2.txt"))
	system(paste0("cat ", outdir, "/Region_*/top_snp_prob_50.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", outdir, "/FINAL_top_snp_prob_50.txt"))
}

##===================
## function - part of summarization - writes fine-mapped SNPs
##===================

Write_Finemap_SNP_Info <- function(BASEDIR, currloci, lociindex, snpdata, CurrRegionOutDir, NUMCAUSALSNP, datatype, finemapmethod) {

	## the snpdata has "beta" in 8th field and "se" in 9'th field
	## estimate p-value
	snpdata$pval <- pnorm((snpdata$beta / snpdata$se), lower.tail=F)
	## sort the snpdata based on decreasing posterior probability (11th field)
	snpdata <- snpdata[order(-snpdata$prob), ]

	## file to contain SNPs sorted by probability, 
	## where the cumulative probability just reaches 95%
	outfile_snp_95pct <- paste0(CurrRegionOutDir, '/snp_cum_prob_95.txt')
	## file to contain top SNPs having probability >= 0.5
	outfile_top_snps_prob_50pct <- paste0(CurrRegionOutDir, '/top_snp_prob_50.txt')
	## file to contain top SNPs having log10BF >= 2 (considered to be causal)
	outfile_top_snps_BF_gt2 <- paste0(CurrRegionOutDir, '/top_snp_log10BF_geq_2.txt')

	## file to contain all the top SNPs in the credible sets
	outfile_top_snps_credible_set <- paste0(CurrRegionOutDir, '/top_snp_credible_set.txt')

	## file to summarize the credible sets (top)
	outfile_top_snps_credible_set_summary <- paste0(CurrRegionOutDir, '/summary_credible_set.txt')

	# compute cumulative posterior probabilities (based on the 11th field) 
	## option 1 - normalized with respect to the total posterior prob
	# cum_prob <- cumsum(snpdata[, 11]) / sum(snpdata[, 11])
	cum_prob <- cumsum(snpdata$prob)

	## SNPs having cumulative posterior probabilities up to 95%
	idx_95 <- min(which(cum_prob >= 0.95))
	## top SNPs with individual posterior probability >= 0.5
	top_snp_idx_prob_50 <- which(snpdata$prob >= 0.5)
	## top SNPs with individual log10BF >= 2
	top_snp_idx_logBF_2 <- which(snpdata$log10bf >= 2)

	if (idx_95 > 0) {
		outdf <- cbind.data.frame(data.frame(regionID=rep(lociindex, idx_95), Loci=rep(currloci, idx_95)), snpdata[1:idx_95, ])
		write.table(outdf, outfile_snp_95pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
	if (length(top_snp_idx_prob_50) > 0) {
		outdf <- cbind.data.frame(data.frame(regionID=rep(lociindex, length(top_snp_idx_prob_50)), Loci=rep(currloci, length(top_snp_idx_prob_50))), snpdata[top_snp_idx_prob_50, ])
		write.table(outdf, outfile_top_snps_prob_50pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
	if (length(top_snp_idx_logBF_2) > 0) {
		outdf <- cbind.data.frame(data.frame(regionID=rep(lociindex, length(top_snp_idx_logBF_2)), Loci=rep(currloci, length(top_snp_idx_logBF_2))), snpdata[top_snp_idx_logBF_2, ])
		write.table(outdf, outfile_top_snps_BF_gt2, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
	cat(sprintf("\n *** Data type : %s Analyzing %s output ---  Total SNPs : %s Number of SNPs with cumulative prob 0.95 : %s Number of top SNPs with posterior probability >= 0.5 : %s Number of top SNPs with log10BF >= 2 : %s ", datatype, finemapmethod, nrow(snpdata), idx_95, length(top_snp_idx_prob_50), length(top_snp_idx_logBF_2)))

	## output files from FINEMAP (.cred)
	## space-delimited text file. It contains the 95% credible sets for each causal signal in the genomic region. 
	## j: number of causal SNPs (and indicates corresponding cred file)			
	rsIDVec <- c()
	for (j in 1:NUMCAUSALSNP) {
		## j: number of causal SNPs for this credible set
		if (datatype == "GWAS") {
			credfile <- paste0(BASEDIR, '/Region', lociindex, '.cred', j)
		} else {
			credfile <- paste0(BASEDIR, '/', currloci, '/Region', lociindex, '.cred', j)
		}
		cat(sprintf("\n ********* %s --- processing credfile : %s ", finemapmethod, credfile))
		if (file.exists(credfile)) {
			creddata <- read.table(credfile, header=T, sep=" ", stringsAsFactors=F)
			## select the first row
			## columns 2, 4, ... contain the SNP IDs in this credible set
			CausalSNPVec <- as.vector(as.character(creddata[1, c(seq(2, ncol(creddata), by=2))]))
			cat(sprintf("\n ********* CausalSNPVec : %s ", paste(CausalSNPVec, collapse=",")))
			currDF <- data.frame(regionID=lociindex, Loci=currloci, numCausalSNP=j, CausalSNPList=paste(CausalSNPVec, collapse=","))
			if (length(rsIDVec) == 0) {
				finalDF <- currDF
			} else {
				finalDF <- rbind.data.frame(finalDF, currDF)
			}
			## add the SNP IDs
			rsIDVec <- c(rsIDVec, CausalSNPVec)
		}
	}
	## unique SNPs from all the causal sets
	rsIDVec <- unique(rsIDVec)
	rsID_DF <- data.frame(rsid=rsIDVec)

	mergeDF <- dplyr::inner_join(snpdata, rsID_DF)
	mergeDF <- cbind.data.frame(data.frame(regionID=rep(lociindex, nrow(mergeDF)), Loci=rep(currloci, nrow(mergeDF))), mergeDF)

	## dump the output files
	write.table(mergeDF, outfile_top_snps_credible_set, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	write.table(finalDF, outfile_top_snps_credible_set_summary, row.names=F, col.names=T, sep="\t", quote=F, append=F)

}

##===================
## this function summarizes finemap output
##===================
Summary_Finemap <- function(inpdir, regiondata, inpdir_cond, inpdir_sss, numSNP, type) {

	## output directory for the summary statistics
	SummaryOutDir <- paste0(inpdir, '/Summary')
	system(paste("mkdir -p", SummaryOutDir))

	## directories to contain summary of finemap output, for individual models (cond and sss) 
	SummaryDir_cond <- paste0(SummaryOutDir, '/cond')
	system(paste("mkdir -p", SummaryDir_cond))
	SummaryDir_sss <- paste0(SummaryOutDir, '/sss')
	system(paste("mkdir -p", SummaryDir_sss))

	bool_cond <- FALSE
	bool_sss <- FALSE

	##======= loop to process regions
	for (i in 1:nrow(regiondata)) {
		if (type == "GWAS") {
			currloci <- paste0(regiondata[i, 1], ":", regiondata[i, 2], "-", regiondata[i, 3])
			cat(sprintf("\n\n ==>>> processing GWAS region number : %s loci : %s ", i, currloci))
		} else {
			currloci = regiondata[i, 1]
			cat(sprintf("\n\n ==>>> processing eGene : %s ", currloci))
		}

		##==========
		## process the output for "cond"
		##==========

		## output file from FINEMAP (.snp)
		## space-delimited text file. 
		## contains GWAS summary statistics and model-averaged posterior summaries for each SNP one per line.
		if (type == "GWAS") {
			snpfile <- paste0(inpdir_cond, '/Region', i, '.snp')		
		} else {
			snpfile <- paste0(inpdir_cond, '/', currloci, '/Region', i, '.snp')
		}
		if ((file.exists(snpfile)) & (file.info(snpfile)$size > 0)) {
			cat(sprintf("\n\n ********* cond --- processing SNP file : %s ", snpfile))
			snpdata <- read.table(snpfile, header=T, sep=" ", stringsAsFactors=F)			
			if (nrow(snpdata) > 0) {
				CurrRegionOutDir <- paste0(SummaryDir_cond, '/', currloci)
				system(paste("mkdir -p", CurrRegionOutDir))
				Write_Finemap_SNP_Info(inpdir_cond, currloci, i, snpdata, CurrRegionOutDir, numSNP, type, "cond")
				bool_cond <- TRUE
			}
		}

		##==========
		## process the output for "sss"
		##==========
		## output file from FINEMAP (.snp)
		## space-delimited text file. 
		## contains the GWAS summary statistics and model-averaged posterior summaries for each SNP one per line.
		if (type == "GWAS") {
			snpfile <- paste0(inpdir_sss, '/Region', i, '.snp')		
		} else {
			snpfile <- paste0(inpdir_sss, '/', currloci, '/Region', i, '.snp')		
		}
		if ((file.exists(snpfile)) & (file.info(snpfile)$size > 0)) {
			cat(sprintf("\n\n ********* sss --- processing SNP file : %s ", snpfile))
			snpdata <- read.table(snpfile, header=T, sep=" ", stringsAsFactors=F)			
			if (nrow(snpdata) > 0) {
				CurrRegionOutDir <- paste0(SummaryDir_sss, '/', currloci)				
				system(paste("mkdir -p", CurrRegionOutDir))
				Write_Finemap_SNP_Info(inpdir_sss, currloci, i, snpdata, CurrRegionOutDir, numSNP, type, "sss")
				bool_sss <- TRUE						
			}
		}	
	}

	## consolidate all the outputs
	if (bool_cond == TRUE) {
		Consolidate_Summary_Output(SummaryDir_cond)
	}

	if (bool_sss == TRUE) {
		Consolidate_Summary_Output(SummaryDir_sss)	
	}
}

##=============
## function to parse GWAS data
##=============
Parse_GWAS <- function(config) {
	GWASData <- data.table::fread(config$Stat$Filename, header=T)

	collist_gwas <- c(as.integer(config$Stat$chrCol), as.integer(config$Stat$posCol))	#, as.integer(config$Stat$allele1Col), as.integer(config$Stat$allele2Col))
	colnameslist_gwas <- c('chr', 'pos')	#, 'allele1', 'allele2')
	if (length(config$Stat$rsIDCol) > 0) {
		collist_gwas <- c(collist_gwas, as.integer(config$Stat$rsIDCol))
		colnameslist_gwas <- c(colnameslist_gwas, 'rsID')
	}
	collist_gwas <- c(collist_gwas, as.integer(config$Stat$pValCol))
	colnameslist_gwas <- c(colnameslist_gwas, 'pval_GWAS')	
	if (length(config$Stat$betaCol) == 0) {
		GWASData$beta <- log(as.numeric(GWASData[, as.integer(config$Stat$ORCol)]))
		collist_gwas <- c(collist_gwas, ncol(GWASData))
	} else {
		collist_gwas <- c(collist_gwas, as.integer(config$Stat$betaCol))
	}
	colnameslist_gwas <- c(colnameslist_gwas, 'beta_GWAS')
	if (length(config$Stat$SECol) > 0) {
		collist_gwas <- c(collist_gwas, as.integer(config$Stat$SECol))
		colnameslist_gwas <- c(colnameslist_gwas, 'SE_gwas')
	}
	if (length(config$Stat$AFCol) > 0) {
		collist_gwas <- c(collist_gwas, as.integer(config$Stat$AFCol))
		colnameslist_gwas <- c(colnameslist_gwas, 'AF_GWAS')
	}
	if (length(config$Stat$SampleSizeCol) > 0) {
		collist_gwas <- c(collist_gwas, as.integer(config$Stat$SampleSizeCol))
		colnameslist_gwas <- c(colnameslist_gwas, 'N_GWAS')
	}

	cat(sprintf("\n collist_gwas : %s ", paste(collist_gwas, collapse=" ")))
	cat(sprintf("\n colnameslist_gwas : %s ", paste(colnameslist_gwas, collapse=" ")))

	GWASData <- GWASData[, c(collist_gwas)]
	colnames(GWASData) <- colnameslist_gwas

	## check data types
	GWASData$pos <- as.integer(GWASData$pos)
	GWASData$beta_GWAS <- as.numeric(GWASData$beta_GWAS)
	if (length(config$Stat$SECol) > 0) {
		GWASData$SE_gwas <- as.numeric(GWASData$SE_gwas)
	}
	GWASData$pval_GWAS <- as.numeric(GWASData$pval_GWAS)	
	if (length(config$Stat$AFCol) > 0) {
		GWASData$AF_GWAS <- as.numeric(GWASData$AF_GWAS)
	}	
	if (length(config$Stat$SampleSizeCol) > 0) {
		GWASData$N_GWAS <- as.integer(GWASData$N_GWAS)
	}

	cat(sprintf("\n\n read GWASData - number of entries : %s ", nrow(GWASData)))

	##====== filter rows having any NA entries
	GWASData <- GWASData[complete.cases(GWASData), ]
	cat(sprintf("\n\n after filtering out NA entries - number of entries : %s ", nrow(GWASData)))

	## if the standard error in GWAS data is not provided, we can compute them
	## using p-value (2-tailed distribution) and beta
	## check: https://www.biostars.org/p/431875/
	if (length(config$Stat$SECol) == 0) {
		GWASData$SE_gwas <- abs(GWASData$beta_GWAS / qnorm(GWASData$pval_GWAS / 2))
	}

	## check the "SE_gwas" column - there should not be any entry with 0 standard error
	## use a small floating point value (say 0.00001)
	idx <- which(GWASData$SE_gwas == 0)
	if (length(idx) > 0) {
		GWASData[idx, c("SE_gwas")] <- DUMMY_SE_VAL
	}
	## if the GWAS data has chromosomes with numbers then append the "chr"
	if (grepl("chr", GWASData[1,1]) == FALSE) {
		GWASData[,1] <- paste0("chr", GWASData[,1])
	}
	## check the "pval_GWAS" column - there should not be any entry with 0
	## use a small floating point value 
	idx <- which(GWASData$pval_GWAS == 0)
	if (length(idx) > 0) {
		GWASData[idx, c("pval_GWAS")] <- DUMMY_P_VAL
	}	

	cat(sprintf("\n *** Number of reference GWAS SNPs: %s ", nrow(GWASData)))
	return(GWASData)
}

##=============
## function to parse eQTL data
##=============
Parse_eQTL <- function(config) {

	eQTLData <- data.table::fread(config$Stat$Filename, header=T)
	collist_eqtl <- c(as.integer(config$Stat$chrCol), as.integer(config$Stat$posCol))
	colnameslist_eqtl <- c('chr', 'pos')
	if (length(config$Stat$geneIDCol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$geneIDCol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'GeneID')
	}
	if (length(config$Stat$geneNameCol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$geneNameCol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'GeneName')
	}
	if (length(config$Stat$rsIDCol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$rsIDCol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'rsID')	
	}
	collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$betaCol))
	colnameslist_eqtl <- c(colnameslist_eqtl, 'beta_eQTL')	
	if (length(config$Stat$SECol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$SECol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'SE_eQTL')
	}
	if (length(config$Stat$AFCol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$AFCol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'AF_eQTL')
	}
	collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$pValCol))
	colnameslist_eqtl <- c(colnameslist_eqtl, 'pval_eQTL')
	# if (length(config$Stat$FDRCol) > 0) {
	# 	collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$FDRCol))
	# 	colnameslist_eqtl <- c(colnameslist_eqtl, 'FDR_eQTL')
	# }
	if (length(config$Stat$SampleSizeCol) > 0) {
		collist_eqtl <- c(collist_eqtl, as.integer(config$Stat$SampleSizeCol))
		colnameslist_eqtl <- c(colnameslist_eqtl, 'N_eQTL')
	}	

	# ## read if there are other columns specified
	# if (length(config$OtherCol) > 0) {
	# 	OtherColList <- as.integer(unlist(strsplit(config$OtherCol, "[,:]")))
	# 	OtherColHeaderList <- unlist(strsplit(config$OtherColHeader, "[,:]"))
	# 	collist_eqtl <- c(collist_eqtl, OtherColList)
	# 	colnameslist_eqtl <- c(colnameslist_eqtl, OtherColHeaderList)
	# }	
	
	eQTLData <- eQTLData[, c(collist_eqtl)]
	colnames(eQTLData) <- colnameslist_eqtl	

	## check data types
	eQTLData$pos <- as.integer(eQTLData$pos)
	eQTLData$beta_eQTL <- as.numeric(eQTLData$beta_eQTL)
	if (length(config$Stat$SECol) > 0) {
		eQTLData$SE_eQTL <- as.numeric(eQTLData$SE_eQTL)
	}
	if (length(config$Stat$AFCol) > 0) {
		eQTLData$AF_eQTL <- as.numeric(eQTLData$AF_eQTL)
	}
	eQTLData$pval_eQTL <- as.numeric(eQTLData$pval_eQTL)
	# if (length(config$Stat$FDRCol) > 0) {
	# 	eQTLData$FDR_eQTL <- as.numeric(eQTLData$FDR_eQTL)
	# }
	if (length(config$Stat$SampleSizeCol) > 0) {
		eQTLData$N_eQTL <- as.integer(eQTLData$N_eQTL)
	}

	cat(sprintf("\n\n read eQTLData - number of entries : %s ", nrow(eQTLData)))

	##====== filter rows having any NA entries
	eQTLData <- eQTLData[complete.cases(eQTLData), ]
	cat(sprintf("\n\n after filtering out NA entries - number of entries : %s ", nrow(eQTLData)))

	## check for any zero entries in "pval_eQTL" column
	## and replace them with dummy floating point values
	idx <- which(eQTLData$pval_eQTL == 0)
	if (length(idx) > 0) {
		eQTLData[idx, c("pval_eQTL")] <- DUMMY_P_VAL
	}	

	## if the standard error in eQTL data is not provided, we can compute them
	## using p-value (2-tailed distribution) and beta
	## check: https://www.biostars.org/p/431875/
	if (length(config$Stat$SECol) == 0) {
		eQTLData$SE_eQTL <- abs(eQTLData$beta_eQTL / qnorm(eQTLData$pval_eQTL / 2))
	}

	## check for any zero entries in "SE_eQTL" column
	## and replace them with dummy floating point values
	idx <- which(eQTLData$SE_eQTL == 0)
	if (length(idx) > 0) {
		eQTLData[idx, c("SE_eQTL")] <- DUMMY_SE_VAL
	}

	## if the eQTL data has chromosomes with numbers then append the "chr"
	if (grepl("chr", eQTLData[1,1]) == FALSE) {
		eQTLData[,1] <- paste0("chr", eQTLData[,1])
	}
	## remove substrings from dots after gene ID
	eQTLData[,3] <- gsub("\\..*","",eQTLData[,3])

	cat(sprintf("\n *** Number of reference eQTLs: %s ", nrow(eQTLData)))
	return(eQTLData)
}


# *************************************************************
# ******** input parameters ***********
# *************************************************************
args <- commandArgs(TRUE)
configfile <- args[1]

## read the configuration file
config <- yaml.load_file(configfile)

OutDir <- config$Path$OutDir
system(paste("mkdir -p", OutDir))

NUMCAUSALSNP <- as.integer(config$Finemap$NUMCAUSALSNP)
NUMTHREAD <- as.integer(config$Finemap$NUMTHREAD)
ldstoreexec <- config$Path$ldstoreexec
finemapexec <- config$Path$finemapexec
GENOTYPEDIR <- config$Path$GENOTYPEDIR
samplecount <- as.integer(config$Path$samplecount)

## output log file
logfile <- paste0(OutDir, '/Out_Summary_', tools::file_path_sans_ext(basename(config$Stat$Filename)), '.log')
sink(logfile)

## folder: inputs formatted for FINEMAP 
FINEMAPInpDir <- paste0(OutDir, '/FINEMAP_INPUT')
system(paste("mkdir -p", FINEMAPInpDir))

## output directory for FINEMAP when --cond (stepwise conditioning) option is used
BaseOutDir_FINEMAP_cond <- paste0(OutDir, '/FINEMAP_OUTPUT/cond')
system(paste("mkdir -p", BaseOutDir_FINEMAP_cond))

## output directory for FINEMAP when --sss (shotgun stochastic search) option is used
BaseOutDir_FINEMAP_sss <- paste0(OutDir, '/FINEMAP_OUTPUT/sss')
system(paste("mkdir -p", BaseOutDir_FINEMAP_sss))

if (config$Finemap$Study == "eQTL") {
	eQTLData <- Parse_eQTL(config)

	## Step 1: run LDStore and FINEMAP
	Run_Finemap_eQTL(eQTLData, OutDir, FINEMAPInpDir, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, samplecount, NUMTHREAD, GENOTYPEDIR, ldstoreexec, finemapexec)

	## Step 2: Summarize finemap output
	if ("GeneName" %in% colnames(eQTLData)) {
		eQTLData_Gene <- data.frame(Gene=unique(eQTLData$GeneName))
	} else {
		eQTLData_Gene <- data.frame(Gene=unique(eQTLData$GeneID))
	}
	Summary_Finemap(OutDir, eQTLData_Gene, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, 'eQTL')

} else if (config$Finemap$Study == "GWAS") {
	GWASData <- Parse_GWAS(config)
	OFFSET <- as.integer(config$Finemap$OFFSET)
	
	## Step 1: extracting GWAS regions
	Input_GWAS_Region_File <- paste0(OutDir, '/GWAS_Regions.txt')	
	Extract_GWAS_Regions(GWASData, OFFSET, Input_GWAS_Region_File)

	## Step 2: run LDStore and FINEMAP
	Run_Finemap_GWAS(FINEMAPInpDir, OutDir, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, GENOTYPEDIR, ldstoreexec, finemapexec, samplecount, NUMCAUSALSNP, NUMTHREAD)

	## Step 3: Summarize finemap output
	regiondata <- read.table(Input_GWAS_Region_File, header=F, sep="\t", stringsAsFactors=F)
	cat(sprintf("\n Number of GWAS regions: %s ", nrow(regiondata)))	
	Summary_Finemap(OutDir, regiondata, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, "GWAS")

}

cat(sprintf("\n\n\n **** Session information **** \n\n\n"))
sessionInfo()
