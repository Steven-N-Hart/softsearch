#!/usr/bin/perl
open(VCF,"$ARGV[0]")||die "Usage: <VCF> <Annotation.bed>\n\n\t\t The annotation BED should be of exons\n";

$bedtools=`which intersectBed`;
if(!$bedtools){die "Requires Bedtools in path\n\n"}
if(!$ARGV[1]){die "Usage: <VCF> <Annotation.bed>\n\n";}

while (<VCF>){
	if($_=~/^#/){print;next}
	chomp;
	@data=split(/\t/,$_);
	#Get left pair information
	$chr1=$data[0];
        $pos1a=$data[1]-1;
        $pos1b=$data[1];
	#Get right pair information
	$data[4]=~s/[ACTGactghr\[\]]//g;#$data[4]=~s/hr/chr/;
	@pos2=split(/:/,$data[4]);
	$chr2="chr";
	$chr2.=$pos2[0];
	$pos2a=$pos2[1]-1;
	$pos2b=$pos2[1];
	#Now get left side annotations
	#
#	print "LEFT=get_anno($chr1,$pos1a,$pos1b)\n";
	$left_gene=get_anno($chr1,$pos1a,$pos1b);
        #print "RIGHT=get_anno($chr2,$pos2a,$pos2b)\n";
        $right_gene=get_anno($chr2,$pos2a,$pos2b);
	print "$_\t$left_gene\t$right_gene\n";
}

close VCF;

sub get_anno(){
	my ($chr,$pos1,$pos2)=@_;
	my $cmd="perl -e 'if(($chr)&&($pos1)&&($pos2)){print join(\"\\t\",$chr,$pos1,$pos2).\"\\n\"}'|intersectBed -a $ARGV[1] -b stdin|cut -f4|head -1";
#print $cmd."\n";print "###########\n";
	$result=`$cmd`;
	$result=~chomp;$result=~s/\n//;
	if(!$result){$result="NA"};
	return $result;
}
