#!/bin/sh
#$ -V
#$ -cwd
#$ -q 1-day
#$ -m ae
#$ -M hart.steven@mayo.edu
#$ -l h_vmem=1G
#$ -l h_stack=10M
BAM=$1
TARGET_BED=$2
SAMPLE_NUMBER=$3

#cat $HEADER > out.${SAMPLE_NUMBER}.sam
samtools view -L $TARGET_BED $BAM|
 perl -ane '
 next if ($F[10]=~/#/);
 $minSize=1000;
 if( $F[1] & 8 || $F[1] & 4 ||  $F[8] == 0 || abs($F[8]) > $minSize || $F[5] =~/S/){
 $rName=join("","@",@F[0]);
  print join ("\n",$rName,$F[9],"+",@F[10])."\n";
};
 ' >> out.${SAMPLE_NUMBER}.fq
echo "Done with $BAM `date`"
