#!/usr/bin/env Rscript

##====================
## R script to summarize FINEMAP outputs
##====================

library(gridExtra)

## Locuszoom option - reference genome + population + build
RefGenomeStr <- '--pop EUR --build hg19 --source 1000G_March2012'

## Locuszoom option - including reference GWAS catalog
GWASCatalogStr <- '--gwas-cat whole-cat_significant-only'

# Print_Credible_Set <- function(snpdata, regionID, idx, currgwasloci, outfile_credible_snp) {
# 	## first print the credible set SNPs
# 	outdf <- cbind.data.frame(data.frame(regionID=rep(regionID, idx), GWASLoci=rep(currgwasloci, idx)), snpdata[1:idx, ])
# 	write.table(outdf, outfile_credible_snp, row.names=F, col.names=T, sep="\t", quot=F, append=F)
# 	## then construct a one line summary
# 	if (idx > 1) {
# 		currDF <- data.frame(regionID=regionID, GWASLoci=currgwasloci, numCausal=idx, LeadSNP=snpdata[1, 2], PP_LeadSNP=snpdata[1, 11], logBF_LeadSNP=snpdata[1, 12], LinkedSNPs=paste(as.vector(snpdata[2:idx, 2]), collapse=","), PP_LinkedSNPs=paste(as.vector(snpdata[2:idx, 11]), collapse=","), logBF_LinkedSNPs=paste(as.vector(snpdata[2:idx, 12]), collapse=","))
# 	} else {
# 		currDF <- data.frame(regionID=regionID, GWASLoci=currgwasloci, numCausal=idx, LeadSNP=snpdata[1, 2], PP_LeadSNP=snpdata[1, 11], logBF_LeadSNP=snpdata[1, 12], LinkedSNPs='-', PP_LinkedSNPs='-', logBF_LinkedSNPs='-')
# 	}
# 	return(currDF)
# }

Plot_LocusZoom <- function(LocusZoomExec, inpfile, tempfile, prefixstr, inpchr, start, end) {
	## create locuszoom input file
	system(paste0("cut -f4,19 ", inpfile, " > ", tempfile))
	## reference SNP
	tempdata <- read.table(tempfile, header=T, sep="\t", stringsAsFactors=F)
	refsnp <- tempdata[1,1]
	
	## first create a plot with flanking 500 Kb option
	## use the first SNP as the reference SNP
	system(paste0(LocusZoomExec, " --metal ", tempfile, " --markercol 1 --pvalcol 2 --refsnp ", refsnp, " --flank 500kb ", RefGenomeStr, " ", GWASCatalogStr, " --plotonly --no-date --prefix ", prefixstr, "_refSNP"))
	
	## then create a plot with respect to the complete GWAS loci
	system(paste0(LocusZoomExec, " --metal ", tempfile, " --markercol 1 --pvalcol 2 ", RefGenomeStr, " ", GWASCatalogStr, " --chr ", inpchr, " --start ", start, " --end ", end, " --plotonly --no-date --prefix ", prefixstr, "_COMPLETE_LOCI"))
}

##==================
## main code
##==================

args <- commandArgs(TRUE)
## input GWAS regions
Input_GWAS_Region_File <- args[1]
## output directory for FINEMAP when --cond (stepwise conditioning) option is used
BaseOutDir_FINEMAP_cond <- args[2]
## output directory for FINEMAP when --sss (shotgun stochastic search) option is used
BaseOutDir_FINEMAP_sss <- args[3]
## number of causal SNPs tested in the finemap execution
NUMCAUSALSNP <- as.integer(args[4])
## locusZoom executable - if provided
if (length(args) > 4) {
	LocusZoomExec <- args[5]
}

## output directory for the summary statistics
SummaryOutDir <- paste0(dirname(Input_GWAS_Region_File), '/Summary')
system(paste("mkdir -p", SummaryOutDir))

outlogfile <- paste0(SummaryOutDir, '/Finemap_Summary.log')
sink(outlogfile)

## directories to contain summary of finemap output, for individual models (cond and sss) 
Finemap_SummaryDir_cond <- paste0(SummaryOutDir, '/cond')
system(paste("mkdir -p", Finemap_SummaryDir_cond))
Finemap_SummaryDir_sss <- paste0(SummaryOutDir, '/sss')
system(paste("mkdir -p", Finemap_SummaryDir_sss))

regiondata <- read.table(Input_GWAS_Region_File, header=F, sep="\t", stringsAsFactors=F)
cat(sprintf("\n Number of GWAS regions: %s ", nrow(regiondata)))

bool_cond <- FALSE
bool_sss <- FALSE

##======= loop to process regions
for (i in 1:nrow(regiondata)) {
	currgwasloci <- paste0(regiondata[i, 1], ":", regiondata[i, 2], "-", regiondata[i, 3])
	cat(sprintf("\n\n ==>>> processing GWAS region number : %s loci : %s ", i, currgwasloci))

	##==========
	## process the output for "cond"
	##==========

	## output file from FINEMAP (.snp)
	## space-delimited text file. It contains the GWAS summary statistics and model-averaged posterior summaries for each SNP one per line.
	snpfile <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.snp')		
	
	if ((file.exists(snpfile)) & (file.info(snpfile)$size > 0)) {
		cat(sprintf("\n\n ********* cond --- processing SNP file : %s ", snpfile))
		snpdata <- read.table(snpfile, header=T, sep=" ", stringsAsFactors=F)
		
		if (nrow(snpdata) > 0) {
			## the snpdata has "beta" in 8th field and "se" in 9'th field
			## estimate p-value
			snpdata$pval <- pnorm((snpdata[,8] / snpdata[,9]), lower.tail=F)
			## sort the snpdata based on decreasing posterior probability (11th field)
			snpdata <- snpdata[order(-snpdata[,11]), ]

			CurrRegionOutDir <- paste0(Finemap_SummaryDir_cond, '/Region_', i)
			system(paste("mkdir -p", CurrRegionOutDir))

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
			cum_prob <- cumsum(snpdata[, 11])

			## SNPs having cumulative posterior probabilities up to 95%
			idx_95 <- min(which(cum_prob >= 0.95))
			## top SNPs with individual posterior probability >= 0.5
			top_snp_idx_prob_50 <- which(snpdata[, 11] >= 0.5)
			## top SNPs with individual log10BF >= 2
			top_snp_idx_logBF_2 <- which(snpdata[, 12] >= 2)

			if (idx_95 > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, idx_95), GWASLoci=rep(currgwasloci, idx_95)), snpdata[1:idx_95, ])
				write.table(outdf, outfile_snp_95pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			if (length(top_snp_idx_prob_50) > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, length(top_snp_idx_prob_50)), GWASLoci=rep(currgwasloci, length(top_snp_idx_prob_50))), snpdata[top_snp_idx_prob_50, ])
				write.table(outdf, outfile_top_snps_prob_50pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			if (length(top_snp_idx_logBF_2) > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, length(top_snp_idx_logBF_2)), GWASLoci=rep(currgwasloci, length(top_snp_idx_logBF_2))), snpdata[top_snp_idx_logBF_2, ])
				write.table(outdf, outfile_top_snps_BF_gt2, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			cat(sprintf("\n *** Analyzing \"cond\" output ---  Total SNPs : %s Number of SNPs with cumulative prob 0.95 : %s Number of top SNPs with posterior probability >= 0.5 : %s Number of top SNPs with log10BF >= 2 : %s ", nrow(snpdata), idx_95, length(top_snp_idx_prob_50), length(top_snp_idx_logBF_2)))

			## output files from FINEMAP (.cred)
			## space-delimited text file. It contains the 95% credible sets for each causal signal in the genomic region. 
			## j: number of causal SNPs (and indicates corresponding cred file)			
			rsIDVec <- c()
			for (j in 1:NUMCAUSALSNP) {
				## j: number of causal SNPs for this credible set
				credfile <- paste0(BaseOutDir_FINEMAP_cond, '/Region', i, '.cred', j)
				cat(sprintf("\n ********* cond --- processing credfile : %s ", credfile))
				if (file.exists(credfile)) {
					creddata <- read.table(credfile, header=T, sep=" ", stringsAsFactors=F)
					## select the first row
					## columns 2, 4, ... contain the SNP IDs in this credible set
					CausalSNPVec <- as.vector(as.character(creddata[1, c(seq(2, ncol(creddata), by=2))]))
					cat(sprintf("\n ********* CausalSNPVec : %s ", paste(CausalSNPVec, collapse=",")))
					currDF <- data.frame(regionID=i, GWASLoci=currgwasloci, numCausalSNP=j, CausalSNPList=paste(CausalSNPVec, collapse=","))
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
			mergeDF <- cbind.data.frame(data.frame(regionID=rep(i, nrow(mergeDF)), GWASLoci=rep(currgwasloci, nrow(mergeDF))), mergeDF)

			## dump the output files
			write.table(mergeDF, outfile_top_snps_credible_set, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			write.table(finalDF, outfile_top_snps_credible_set_summary, row.names=F, col.names=T, sep="\t", quote=F, append=F)			

			##=============
			## now plot the SNPs in locusZoom
			##=============
			if (length(args) > 4) {
				temp_locuszoom_inp_file <- paste0(CurrRegionOutDir, '/temp_locuszoom_input.txt')			
				inpchr <- snpdata[1,5]
				if (grepl('chr', inpchr) == FALSE) {
					inpchr <- paste0('chr', inpchr)
				}
				if (file.exists(outfile_snp_95pct)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/snp_cum_prob_95')
					Plot_LocusZoom(LocusZoomExec, outfile_snp_95pct, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_prob_50pct)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snp_prob_50')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_prob_50pct, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_BF_gt2)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snps_BF_gt2')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_BF_gt2, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_credible_set)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snps_credible_set')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_credible_set, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(temp_locuszoom_inp_file)) {
					system(paste("rm", temp_locuszoom_inp_file))
				}
			}
			##=======
			## end plot
			##=======

			## boolean variable set
			bool_cond <- TRUE
		}
	}	

	##==========
	## process the output for "sss"
	##==========
	## output file from FINEMAP (.snp)
	## space-delimited text file. It contains the GWAS summary statistics and model-averaged posterior summaries for each SNP one per line.
	snpfile <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.snp')		

	if ((file.exists(snpfile)) & (file.info(snpfile)$size > 0)) {
		cat(sprintf("\n\n ********* sss --- processing SNP file : %s ", snpfile))
		snpdata <- read.table(snpfile, header=T, sep=" ", stringsAsFactors=F)
		
		if (nrow(snpdata) > 0) {
			## the snpdata has "beta" in 8th field and "se" in 9'th field
			## estimate p-value
			snpdata$pval <- pnorm((snpdata[,8] / snpdata[,9]), lower.tail=F)

			## sort the snpdata based on decreasing posterior probability (11th field)
			snpdata <- snpdata[order(-snpdata[,11]), ]

			CurrRegionOutDir <- paste0(Finemap_SummaryDir_sss, '/Region_', i)
			system(paste("mkdir -p", CurrRegionOutDir))

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

			## compute cumulative posterior probabilities (based on the 11th field)
			## option 1 - normalized with respect to the total posterior prob
			# cum_prob <- cumsum(snpdata[, 11]) / sum(snpdata[, 11])
			cum_prob <- cumsum(snpdata[, 11])

			## SNPs having cumulative posterior probabilities up to 95%
			idx_95 <- min(which(cum_prob >= 0.95))
			## top SNPs with individual posterior probability >= 0.5
			top_snp_idx_prob_50 <- which(snpdata[, 11] >= 0.5)
			## top SNPs with individual log10BF >= 2
			top_snp_idx_logBF_2 <- which(snpdata[, 12] >= 2)

			if (idx_95 > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, idx_95), GWASLoci=rep(currgwasloci, idx_95)), snpdata[1:idx_95, ])
				write.table(outdf, outfile_snp_95pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			if (length(top_snp_idx_prob_50) > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, length(top_snp_idx_prob_50)), GWASLoci=rep(currgwasloci, length(top_snp_idx_prob_50))), snpdata[top_snp_idx_prob_50, ])
				write.table(outdf, outfile_top_snps_prob_50pct, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			if (length(top_snp_idx_logBF_2) > 0) {
				outdf <- cbind.data.frame(data.frame(regionID=rep(i, length(top_snp_idx_logBF_2)), GWASLoci=rep(currgwasloci, length(top_snp_idx_logBF_2))), snpdata[top_snp_idx_logBF_2, ])
				write.table(outdf, outfile_top_snps_BF_gt2, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			}
			cat(sprintf("\n *** Analyzing \"sss\" output ---  Total SNPs : %s Number of SNPs with cumulative prob 0.95 : %s Number of top SNPs with posterior probability >= 0.5 : %s Number of top SNPs with log10BF >= 2 : %s ", nrow(snpdata), idx_95, length(top_snp_idx_prob_50), length(top_snp_idx_logBF_2)))

			## output files from FINEMAP (.cred)
			## space-delimited text file. It contains the 95% credible sets for each causal signal in the genomic region. 
			## j: number of causal SNPs (and indicates corresponding cred file)			
			rsIDVec <- c()
			for (j in 1:NUMCAUSALSNP) {
				## j: number of causal SNPs for this credible set
				credfile <- paste0(BaseOutDir_FINEMAP_sss, '/Region', i, '.cred', j)
				cat(sprintf("\n ********* sss --- processing credfile : %s ", credfile))
				if (file.exists(credfile)) {
					creddata <- read.table(credfile, header=T, sep=" ", stringsAsFactors=F)
					## select the first row
					## columns 2, 4, ... contain the SNP IDs in this credible set
					CausalSNPVec <- as.vector(as.character(creddata[1, c(seq(2, ncol(creddata), by=2))]))
					cat(sprintf("\n ********* CausalSNPVec : %s ", paste(CausalSNPVec, collapse=",")))
					currDF <- data.frame(regionID=i, GWASLoci=currgwasloci, numCausalSNP=j, CausalSNPList=paste(CausalSNPVec, collapse=","))
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
			rsID_DF <- data.frame(rsid=rsIDVec, stringsAsFactors=F)

			mergeDF <- dplyr::inner_join(snpdata, rsID_DF)
			mergeDF <- cbind.data.frame(data.frame(regionID=rep(i, nrow(mergeDF)), GWASLoci=rep(currgwasloci, nrow(mergeDF))), mergeDF)

			## dump the output files
			write.table(mergeDF, outfile_top_snps_credible_set, row.names=F, col.names=T, sep="\t", quote=F, append=F)
			write.table(finalDF, outfile_top_snps_credible_set_summary, row.names=F, col.names=T, sep="\t", quote=F, append=F)			

			##=============
			## now plot the SNPs in locusZoom
			##=============
			if (length(args) > 4) {
				temp_locuszoom_inp_file <- paste0(CurrRegionOutDir, '/temp_locuszoom_input.txt')			
				inpchr <- snpdata[1,5]
				if (grepl('chr', inpchr) == FALSE) {
					inpchr <- paste0('chr', inpchr)
				}
				if (file.exists(outfile_snp_95pct)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/snp_cum_prob_95')
					Plot_LocusZoom(LocusZoomExec, outfile_snp_95pct, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_prob_50pct)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snp_prob_50')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_prob_50pct, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_BF_gt2)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snps_BF_gt2')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_BF_gt2, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(outfile_top_snps_credible_set)) {
					## prefix string for locuszoom plot
					prefixstr <- paste0(CurrRegionOutDir, '/top_snps_credible_set')
					Plot_LocusZoom(LocusZoomExec, outfile_top_snps_credible_set, temp_locuszoom_inp_file, prefixstr, inpchr, regiondata[i,2], regiondata[i,3])
				}
				if (file.exists(temp_locuszoom_inp_file)) {
					system(paste("rm", temp_locuszoom_inp_file))
				}
			}
			##=======
			## end plot
			##=======

			## boolean variable set
			bool_sss <- TRUE						
		}
	}	
}

##=======
## consolidate all the outputs
##=======
if (bool_cond == TRUE) {
	system(paste0("cat ", Finemap_SummaryDir_cond, "/Region_*/snp_cum_prob_95.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_cond, "/FINAL_snp_cum_prob_95.txt"))
	system(paste0("cat ", Finemap_SummaryDir_cond, "/Region_*/summary_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_cond, "/FINAL_summary_credible_set.txt"))
	system(paste0("cat ", Finemap_SummaryDir_cond, "/Region_*/top_snp_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_cond, "/FINAL_top_snp_credible_set.txt"))
	system(paste0("cat ", Finemap_SummaryDir_cond, "/Region_*/top_snp_log10BF_geq_2.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_cond, "/FINAL_top_snp_log10BF_geq_2.txt"))
	system(paste0("cat ", Finemap_SummaryDir_cond, "/Region_*/top_snp_prob_50.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_cond, "/FINAL_top_snp_prob_50.txt"))
}

if (bool_sss == TRUE) {
	system(paste0("cat ", Finemap_SummaryDir_sss, "/Region_*/snp_cum_prob_95.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_sss, "/FINAL_snp_cum_prob_95.txt"))
	system(paste0("cat ", Finemap_SummaryDir_sss, "/Region_*/summary_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_sss, "/FINAL_summary_credible_set.txt"))
	system(paste0("cat ", Finemap_SummaryDir_sss, "/Region_*/top_snp_credible_set.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_sss, "/FINAL_top_snp_credible_set.txt"))
	system(paste0("cat ", Finemap_SummaryDir_sss, "/Region_*/top_snp_log10BF_geq_2.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_sss, "/FINAL_top_snp_log10BF_geq_2.txt"))
	system(paste0("cat ", Finemap_SummaryDir_sss, "/Region_*/top_snp_prob_50.txt | awk \'!((NR>1) && ($1==\"regionID\"))\' - > ", Finemap_SummaryDir_sss, "/FINAL_top_snp_prob_50.txt"))
}

