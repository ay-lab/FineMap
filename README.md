Fine mapping of GWAS and eQTL summary statistics
-------------------------------------------------

Sourya Bhattacharyya

La Jolla Institute for Immunology, La Jolla, CA 92037, USA

----------------------

Implements a wrapper for the FINEMAP package (http://www.christianbenner.com/), to generate the fine mapped variants corresponding to an input GWAS or eQTL summary statistics


Prerequisites
===============

1. 	Download the 1000G Genotype data (we've used hg19 as the reference genome).

	Follow the tutorial: https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html. 

	The path containing this 1000G genotype data needs to be provided as a configuration parameter (check the parameter *GENOTYPEDIR* in the configuration files, described below).

2. 	Install these R libraries: data.table, dplyr, GenomicRanges, gridExtra


3. 	Download the package LDSTORE2 (binary executable) from http://www.christianbenner.com/. 

	The path of this executable needs to be provided as a configuration parameter (check the parameter *ldstoreexec* in the configuration files, described below)

4. 	Download the package FINEMAP (binary executable) from http://www.christianbenner.com/. 

	The path of this executable needs to be provided as a configuration parameter (check the parameter *finemapexec* in the configuration files, described below)


Fine-mapping of eQTL summary statistics
==========================================

First edit the configuration file *configfile_eQTL.yaml*.

A sample eQTL smmary statistics is provided in the file *Data/Example_eQTL_data.txt* for the users.

User needs to edit the configuration file according to the file / directory paths, contents of the summary statistics file.

	.. Note ..

		1. In the configuration file, user needs to provide at least one of *betaCol* (column containing the beta statistics) or *ORCol* (column containing the odds ratio).

		2. Either user needs to provide the *SampleSizeCol* (column containing the sample size information) or put the integer value in the field *TotalSample* (total number of samples / donors used in the eQTL study).

		3. The parameter *geneIDCol* (column containing the gene ID, preferably in Ensemble ID format) is mandatory.


Once the configuration file is edited, user can use the script *finemap_script.sh* to execute the fine-mapping.


Fine-mapping of GWAS summary statistics
==========================================

First edit the configuration file *configfile_GWAS.yaml*.

A sample GWAS smmary statistics is provided in the file *Data/Example_GWAS_data.txt* for the users.

User needs to edit the configuration file according to the file / directory paths, contents of the summary statistics file.

	.. Note ..

		1. In the configuration file, user needs to provide at least one of *betaCol* (column containing the beta statistics) or *ORCol* (column containing the odds ratio).

		2. Either user needs to provide the *SampleSizeCol* (column containing the sample size information) or put the integer value in the field *TotalSample* (total number of samples / donors used in the eQTL study).

		3. The parameters *chrCol* (column containing the chromosome) and *posCol* (column containing the SNP position) are mandatory. 

		4. User should provide the parameter *AFCol* (column containing the allele frequency information), if available.

		5. The *OFFSET* parameter specifies the window surrounding a lead GWAS SNP used to define a GWAS locus.


Once the configuration file is edited, user can use the script *finemap_script.sh* to execute the fine-mapping.



Output
========

With respect to the specified output directory (parameter *OutDir* in either of these configuration files), check the folders *cond* (output of stepwise conditioning) and *sss* (output of shotgun stochastic search).

	1. *FINAL_summary_credible_set.txt* : File containing the GWAS loci, credible causal sets and the corresponding SNPs.

	2. *FINAL_top_snp_credible_set.txt* : All the SNPs listed in the credible causal sets. 


Contact
==========

For any queries, please create an issue and we'll respond. Otherwise, please e-mail

Sourya Bhattacharyya: sourya@lji.org

Ferhat Ay: ferhatay@lji.org

