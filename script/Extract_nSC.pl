#!/usr/bin/perl -w

use Getopt::Long;

#Initialize values
my (@queries,@HEADER,$samples,@HEADER_OUT,$end,$samp);
GetOptions ("query|q=s" => \$queries);
if(!$queries){die "Usage: FORMAT_extract.pl <VCF> -query nSC 
\n\n";}


open (VCF,"$ARGV[0]") or die "Usage: <VCF>";

while (<VCF>) {
        if($_=~/^##/){print;next}
    chomp;
    @line=split(/\t/,$_);
    if($line[0]=~/^#CH/){
        print join ("\t",@line,$queries)."\n";
	next}
 @FORMAT=split(/:/,$line[8]);
 @SAMPLE=split(/:/,$line[9]);
	for($i=0;$i<@FORMAT;$i++){
	if($FORMAT[$i] =~/^$queries$/){print join ("\t",@line,$SAMPLE[$i])."\n";next}
	}
}
close VCF;
