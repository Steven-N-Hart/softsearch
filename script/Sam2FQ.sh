#!/bin/sh
#$ -V
#$ -cwd
#$ -q 1-day
#$ -m ae
#$ -M hart.steven@mayo.edu
#$ -l h_vmem=20G
#$ -l h_stack=10M
SAM=$1
PICARD_PATH=/projects/bsi/bictools/apps/alignment/picard/1.96
#Seperate single and paired reads
#/usr/local/biotools/java/jdk1.7.0_03/bin/java -Xmx16G -jar $PICARD_PATH/FixMateInformation.jar I=${SAM} O=${SAM}.fixed.sam

samtools view -hS -f 1 ${SAM}.fixed.sam > ${SAM}.paired.sam
samtools view -hS -F 1 ${SAM}.fixed.sam > ${SAM}.unpaired.sam
rm ${SAM}
/usr/local/biotools/java/jdk1.7.0_03/bin/java -Xmx16G -jar $PICARD_PATH/SamToFastq.jar I=${SAM}.paired.sam F=${SAM%%.sam}.1.fq INTERLEAVE=true TMP_DIR=$PWD 
/usr/local/biotools/java/jdk1.7.0_03/bin/java -Xmx16G -jar $PICARD_PATH/SamToFastq.jar I=${SAM}.unpaired.sam F=${SAM%%.sam}.2.fq TMP_DIR=$PWD


