#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=4:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

bash Finemap_Script.sh -C configfile


