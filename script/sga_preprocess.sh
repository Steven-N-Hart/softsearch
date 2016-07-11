#!/bin/sh
#$ -V
#$ -cwd
#$ -q 1-day
#$ -m ae
#$ -M hart.steven@mayo.edu
#$ -l h_vmem=1G
#$ -l h_stack=10M

SGA_PATH=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/sga/bin/
READS_FQ=$1
BAM_NUMBER=$2
$SGA_PATH/sga preprocess -o prepro.${$BAM_NUMBER}.out --dust -p 2 $READS_FQ
