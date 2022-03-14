#!/bin/bash

usage(){
cat << EOF

Options:
   	-C  ConfigFile		Name of the configuration file storing the parameters of FitHiChIP.
EOF
}

while getopts "C:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

echo -e "\n\n ================ Parsing input configuration file ================= \n\n"

## default initialization
LocusZoomExec=''

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			# if there are multiple parameter values (separated by # - old values are kept)
			# then the following operation selects the current one
			paramval=$(echo "$paramval" | awk -F['#\t'] '{print $1}' | tr -d '[:space:]');
			echo -e "Content of $param is $paramval"
			if [ $param == "InpGWASFile" ]; then
				InpGWASFile=$paramval
			fi
			if [ $param == "OutDir" ]; then
				BaseOutDir_GWAS=$paramval
			fi
			if [ $param == "chrcol" ]; then
				chrcol=$paramval
			fi
			if [ $param == "poscol" ]; then
				poscol=$paramval
			fi
			if [ $param == "pvalcol" ]; then
				pvalcol=$paramval
			fi
			if [ $param == "afcol" ]; then
				afcol=$paramval
			fi
			if [ $param == "betacol" ]; then
				betacol=$paramval
			fi
			if [ $param == "secol" ]; then
				secol=$paramval
			fi
			if [ $param == "OFFSET" ]; then
				OFFSET=$paramval
			fi
			if [ $param == "ldstoreexec" ]; then
				ldstoreexec=$paramval
			fi
			if [ $param == "finemapexec" ]; then
				finemapexec=$paramval
			fi
			if [ $param == "GENOTYPEDIR" ]; then
				GENOTYPEDIR=$paramval
			fi
			if [ $param == "samplecount" ]; then
				samplecount=$paramval
			fi
			if [ $param == "NUMCAUSALSNP" ]; then
				NUMCAUSALSNP=$paramval
			fi
			if [ $param == "NUMTHREAD" ]; then
				NUMTHREAD=$paramval
			fi
			if [ $param == "LocusZoomExec" ]; then
				LocusZoomExec=$paramval
			fi
		fi
	fi
done < $ConfigFile

## output directory
mkdir -p ${BaseOutDir_GWAS}

## dump the parameters
ParameterFile=${BaseOutDir_GWAS}'/Parameters.txt'
echo '** Input parameters **' > $ParameterFile
echo 'InpGWASFile : '$InpGWASFile >> $ParameterFile
echo 'BaseOutDir_GWAS : '$BaseOutDir_GWAS >> $ParameterFile
echo 'chrcol : '$chrcol >> $ParameterFile
echo 'poscol : '$poscol >> $ParameterFile
echo 'pvalcol : '$pvalcol >> $ParameterFile
echo 'afcol : '$afcol >> $ParameterFile
echo 'betacol : '$betacol >> $ParameterFile
echo 'secol : '$secol >> $ParameterFile
echo 'OFFSET : '$OFFSET >> $ParameterFile
echo 'ldstoreexec : '$ldstoreexec >> $ParameterFile
echo 'finemapexec : '$finemapexec >> $ParameterFile
echo 'GENOTYPEDIR : '$GENOTYPEDIR >> $ParameterFile
echo 'samplecount : '$samplecount >> $ParameterFile
echo 'NUMCAUSALSNP : '$NUMCAUSALSNP >> $ParameterFile
echo 'NUMTHREAD : '$NUMTHREAD >> $ParameterFile
echo 'LocusZoomExec : '$LocusZoomExec >> $ParameterFile

## identify the directory containing this script
currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir

GWASRegionCodeExec=$currscriptdir'/Extract_GWAS_Regions.R'
LDStoreTestCodeExec=$currscriptdir'/Check_SNPs_Input_Ref.R'
SummaryCodeExec=$currscriptdir'/Summary_FINEMAP_Output.R'

## this folder will contain the inputs formatted for FINEMAP 
FINEMAPInpDir=$BaseOutDir_GWAS'/FINEMAP_INPUT'
mkdir -p ${FINEMAPInpDir}

## output directory for FINEMAP when --cond (stepwise conditioning) option is used
BaseOutDir_FINEMAP_cond=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/cond'
mkdir -p $BaseOutDir_FINEMAP_cond

## output directory for FINEMAP when --sss (shotgun stochastic search) option is used
BaseOutDir_FINEMAP_sss=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/sss'
mkdir -p $BaseOutDir_FINEMAP_sss

## input GWAS regions
Input_GWAS_Region_File=${BaseOutDir_GWAS}'/GWAS_Regions.txt'

## output log file
logfile=${BaseOutDir_GWAS}'/out.log'
echo 'output log' > $logfile

##============
## Step 1: extracting GWAS regions
## use the $OFFSET values
##============
Rscript $GWASRegionCodeExec $InpGWASFile ${BaseOutDir_GWAS} $chrcol $poscol $pvalcol $OFFSET

##============
## Step 2: run LDStore
##============

## create a temporary file
tempfile=${FINEMAPInpDir}'/t1.txt'

##=== process individual GWAS regions
i=0
while [ 1 == 1 ]; do
	i=`expr $i + 1`
	gwasfile=${BaseOutDir_GWAS}'/GWAS_Regions/Region_'${i}'.txt'
	if [[ ! -f $gwasfile ]]; then
		break
	fi
	echo 'processing gwasfile : '$gwasfile >> $logfile
	finemap_z_file=${FINEMAPInpDir}'/Region'${i}'.z'

	## master file which will be used as the input of FINEMAP
	Masterfile=${FINEMAPInpDir}'/input.Region_'${i}'.master'
	echo "z;bgen;bgi;bcor;ld;n_samples;bdose" > ${Masterfile}

	## get the chromosome
	awk '{if (NR==2) {print $0}}' ${gwasfile} > $tempfile	
	currchr=`awk '{if (NR==1) {print $1}}' ${tempfile}`
	echo 'current chromosome : '$currchr >> $logfile
	## get the chromosome number
	sed -i 's/chr//g' ${tempfile}
	currchrnum=`awk '{if (NR==1) {print $1}}' ${tempfile}`
	echo 'current chromosome number : '$currchrnum >> $logfile

	## reference 1000G genotypes - variants are referred by chrNum:pos
	## we are using all the reference panels
	bgenfile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bgen'
	bgenbgifile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bgen.bgi'
	bimfile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bim'

	## use chromosome:position format for variant ID
	## this file will also be useful for FINEMAP
	## note: 6th field requires MAF - so use < 0.5 and convert the allele frequency if required
	## space separated file
	echo "rsid chromosome position maf beta se" > ${finemap_z_file}
	awk -v n="${currchrnum}" -v C="$chrcol" -v S="$poscol" -v F="$afcol" -v B="$betacol" -v E="$secol" '{if (NR>1) {if ($F>0.5) {f=(1-$F)} else {f=$F}; print n":"$S" "n" "$S" "f" "$B" "$E}}'  $gwasfile >> ${finemap_z_file}	

	##========= before running LDStore, check if these SNPs are present in the reference genotype (PLINK) output
	##========= otherwise filter those entries
	n=`cat ${finemap_z_file} | wc -l`
	echo 'before filtering SNPs from FINEMAP input - number of entries : '$n >> $logfile
	Rscript ${LDStoreTestCodeExec} ${finemap_z_file} ${bimfile}
	n=`cat ${finemap_z_file} | wc -l`
	echo 'after filtering SNPs from FINEMAP input - number of entries : '$n >> $logfile

	##======== another filtering - sourya
	##======== sometimes same SNPs (same rsID, position) exist multiple times, even in the original GWAS data - filter them
	awk -F" " '!seen[$1]++' ${finemap_z_file} > $tempfile
	mv $tempfile ${finemap_z_file}
	n=`cat ${finemap_z_file} | wc -l`
	echo 'after removing any duplicate SNPs from FINEMAP input - number of entries : '$n >> $logfile

	##========= write the master file entries
	outbcorfile=${FINEMAPInpDir}'/Region'${i}'.bcor'
	outldfile=${FINEMAPInpDir}'/Region'${i}'.ld'
	bdosefile=${FINEMAPInpDir}'/Region'${i}'.bdose'
	echo ${finemap_z_file}';'${bgenfile}';'${bgenbgifile}';'${outbcorfile}';'${outldfile}';'${samplecount}';'${bdosefile} >> ${Masterfile}

	## now execute ldstore using the generared master file
	${ldstoreexec} --in-files ${Masterfile} --read-only-bgen --write-text --write-bdose --n-threads ${NUMTHREAD} --memory 40

done

##============
## Step 3: run finemap
##============

i=0
while [ 1 == 1 ]; do
	i=`expr $i + 1`
	gwasfile=${BaseOutDir_GWAS}'/GWAS_Regions/Region_'${i}'.txt'
	if [[ ! -f $gwasfile ]]; then
		break
	fi

	## master file containing the finemap execution input
	## when --cond (stepwise conditioning) option is used 
	Masterfile_cond=${BaseOutDir_FINEMAP_cond}'/input_Region_'${i}'.master'
	echo "z;ld;snp;config;cred;log;n_samples" > ${Masterfile_cond}

	## master file containing the finemap execution input
	## when --sss (shotgun stochastic search) option is used 
	Masterfile_sss=${BaseOutDir_FINEMAP_sss}'/input_Region_'${i}'.master'
	echo "z;ld;snp;config;cred;log;n_samples" > ${Masterfile_sss}

	finemap_z_file=${FINEMAPInpDir}'/Region'${i}'.z'
	finemap_ld_file=${FINEMAPInpDir}'/Region'${i}'.ld'

	## construct the master file "Masterfile_cond"
	finemap_out_snp_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.snp'
	finemap_out_config_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.config'
	finemap_out_cred_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.cred'
	finemap_out_log_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.log'
	echo ${finemap_z_file}';'${finemap_ld_file}';'${finemap_out_snp_file}';'${finemap_out_config_file}';'${finemap_out_cred_file}';'${finemap_out_log_file}';'${samplecount} >> ${Masterfile_cond}

	## construct the master file "Masterfile_sss"
	finemap_out_snp_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.snp'
	finemap_out_config_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.config'
	finemap_out_cred_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.cred'
	finemap_out_log_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.log'
	echo ${finemap_z_file}';'${finemap_ld_file}';'${finemap_out_snp_file}';'${finemap_out_config_file}';'${finemap_out_cred_file}';'${finemap_out_log_file}';'${samplecount} >> ${Masterfile_sss}	

	## fine mapping when --cond (stepwise conditioning) option is used 
	${finemapexec} --cond --in-files ${Masterfile_cond} --log --n-causal-snps ${NUMCAUSALSNP}

	## fine mapping when --sss (shotgun stochastic search) option is used 
	${finemapexec} --sss --in-files ${Masterfile_sss} --log --n-causal-snps ${NUMCAUSALSNP} --n-threads ${NUMTHREAD} 

done

##=====================
## step 4: summarize fine mapping output
##=====================
Rscript $SummaryCodeExec ${Input_GWAS_Region_File} ${BaseOutDir_FINEMAP_cond} ${BaseOutDir_FINEMAP_sss} ${NUMCAUSALSNP}


##===========
## now go back to the original working directory
##===========
cd $currworkdir
echo 'Thank you !! Fine mapping is executed.'

