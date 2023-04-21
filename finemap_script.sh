#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=10:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

##====== FINEMAP on sample GWAS data
configfile="configfile_GWAS.yaml"
Rscript Finemap_MAIN.R $configfile

##===== FINEMAP on sample eQTL data
configfile="configfile_eQTL.yaml"
Rscript Finemap_MAIN.R $configfile



