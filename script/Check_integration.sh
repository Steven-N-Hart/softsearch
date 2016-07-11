#!/bin/sh
#$ -V
#$ -cwd
#$ -q 1-day
#$ -m ae
#$ -M hart.steven@mayo.edu
#$ -l h_vmem=8G
#$ -l h_stack=10M
VCF_FILE=$1
x=$2
#VIRAL_SEQDB=/data2/bsi/tertiary/m110344/SoftTile/Mia/BLAST_DB/OBrien/Virus_PCGS.fasta #must me indexed by bwasw
VIRAL_SEQDB=$3
VCF_FILE=$4
#VCF_FILE=final.vcf

set -x 

perl -ane '$dist=100;
$mate=$F[4];
$mate=~s/[A-Z]|\[|\]//g;
@mate=split(/:/,$mate);
$end1a=@F[1]-$dist;
$end1b=@F[1]+$dist;
$end2a=$dist+@mate[1];
$end2b=$dist+@mate[1];
print "@F[0]\t$end1a\t$end1b\n@mate[0]\t$end2a\t$end2b\n"' $VCF_FILE|
sortBed > targets.bed


#100 min
time samtools view -h $x -L targets.bed |awk '(($9==0)&&($11!~/#/)&&($3!~/^chrGL/)&&($3!~/^chrM/))'|perl -ane 'print "\@@F[0]\n@F[9]\n+\n@F[10]\n"' >${x%%.bam}.res.fq
#23 min
time bwa mem -t 4 $VIRAL_SEQDB ${x%%.bam}.res.fq |samtools view -S - |grep gi > ${x%%.bam}.tmp.sam

#find out how many hits there are
cut -f3 ${x%%.bam}.tmp.sam|perl -ne '@_=split(":",$_);@res=split(/_/,@_[1],2);print "@res[1]"' | sort|uniq -c|sort -k1nr|tee ${x%%.bam}.Viral_maps.out |head
#Get the reads mapping to those hits to find out where the integration site is

#Read in the viruses until there is a significant drop off in number of reads (i.e. contributing less than 10%)
perl -ne '@_=split(" ",$_);$i=$_[0]+$i;$j=$_[0];$res=$j/$i;if($res>.1){print "@_[1]\n"}' ${x%%.bam}.Viral_maps.out >${x%%.bam}.to.keep
fgrep -f ${x%%.bam}.to.keep ${x%%.bam}.tmp.sam |cut -f1 >${x%%.bam}.reads

#75min+

time samtools view $x -L targets.bed |
fgrep -f ${x%%.bam}.reads|
awk '{if(($9==0)&&($11!~/#/)&&($3!~/^chrGL/)&&($3!~/^chrM/)&&($3!~/^\*/)){print $3"\t"$4"\t"$4"\t"$1}}'|
tee ${x%%.bam}.unsorted.bed|
sortBed | mergeBed -nms -d 1000|
perl -e 'open (FILE,"$ARGV[0]") or die "cant open file\n\n";
 $SAM="$ARGV[1]";
 $SAM=~chomp;
 while(<FILE>){
chomp;
  @_=split(/\t/,$_);
  @reads=split(/;/,@_[3]);
#print "LINE=$_\nRES=grep $reads[0] $SAM\n";
  $res=`grep $reads[0] $SAM` ;
#  print "AFTER GREP, RES=$res\n";
  if($res){
   @res=split(/\t/,$res);
   print join("\t",@_[0..2],@res[2])."\n"
   }
  };
 close FILE' - ${x%%.bam}.tmp.sam |
perl -ne 's/\|/\t/g;@_=split("\t",$_);print join ("\t",@_[0..2,7])'|
perl -pne 's/_/\t/'|  cut -f4 --complement |
perl -e '
 open (FILE,"$ARGV[0]") or die "cant open file\n\n";
 $SAM="$ARGV[1]";
 while(<FILE>){
  chomp;
  @_=split(/\t/,$_);
  $res=`grep $_[3] $SAM`;
  if($res){
   @res=split(" ",$res);
   $reads[0]=chomp;
   print join("\t",@_[0..4],@res[0])."\n";
  }
 }
close FILE;
' - ${x%%.bam}.Viral_maps.out|
perl -pne 's/_/ /g'> ${x%%.bam}.Virus.integrated.bed

rm ${x%%.bam}.reads ${x%%.bam}.to.keep ${x%%.bam}.tmp.sam ${x%%.bam}.res.fq
echo "DONE with $x"
