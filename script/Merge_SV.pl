#!/usr/bin/perl -w
use Getopt::Long;
use List::Util qw(min max);


#Declare variables
my ($window,$tmpSpace,$usage,$help,$outFile);

GetOptions(
        'v=s{2,}' => \@VCF,
        'o:s' => \$outFile,
        'w:s' => \$window,
		'h|help' => \$help
);

if((!@VCF)||($help)){&usage();exit}


if (!$window) {
    $window=500;
}
if (!$outFile) {
    $outFile="merged.vcf.out";
}
###########################################
# Protect against merging too many results
###########################################
$tmpSpace='temporarySV_merge';
if (-e $tmpSpace) {
    #Delete temp file if it exists
    unlink $tmpSpace;
}
###########################################
#For each VCF, create a BEDPE file
###########################################

open(OUT,">>$tmpSpace") or die "Can't write in this directory\n";
for (my $i=0;$i<@VCF;$i++){
    #print STDERR "opening $VCF[$i]\n";
    open(VCF,$VCF[$i]) or die &usage();
    while (<VCF>) {
        next if ($_=~/^#/);
        chomp;
        @line=split("\t",$_);
        $mate=$line[4];
        $mate=~s/[A-L]|[N-W]|[Z]|\[|\]//g;
        @mate=split(/:/,$mate);
        $end1a=$line[1]-$window;
        $end1b=$line[1]+$window;
        $end2a=$mate[1]-$window;
        $end2b=$mate[1]+$window;
        next if (($end1a<0)||($end2a<0));
        if (($line[0]=~/^chr$/)||($mate[0]=~/^chr$/)) {
            next;
        }
        print OUT "$line[0]\t$end1a\t$end1b\t$mate[0]\t$end2a\t$end2b\n";
        print OUT "$mate[0]\t$end2a\t$end2b\t$line[0]\t$end1a\t$end1b\n";
    }
}
close OUT;

###########################################
#Now merge the BEDPE into a unique BEDPE
###########################################
#Make sure the BEDPE is sorted
#print "Make sure the BEDPE is sorted\n";
my $tmpSpace2=join("",$tmpSpace,".2");
system("cat $tmpSpace|sort -k1,1 -k2,3n -k4,4 -k5,5n -u > $tmpSpace2");
unlink($tmpSpace);

#Create output files for the left and right merged BEDPE
my $tmpSpace3=join("",$tmpSpace,".3");
my $tmpSpace4=join("",$tmpSpace,".4");

open (OUT1,">$tmpSpace3") or die "Cant write in this directory\n";
open (OUT2,">$tmpSpace4") or die "Cant write in this directory\n";

open(BEDPE,"$tmpSpace2") or die "$tmpSpace2 has already been deleted\n";
#Initialize positions
#my ($chr1,$pos2,$pos3,$chr2,$pos3,$pos4);
my (@chr,@pos1,@pos2,@chr2,@pos3,@pos4);
while (<BEDPE>) {
    ($chr1,$pos2,$pos3,$chr2,$pos3,$pos4)=split("\t",$_);
	if(!$Echr1){
	($Echr1,$Epos1,$Epos2,$Echr2,$Epos3,$Epos4)=split("\t",$_);
	}
    while ( 
		 ($chr1 =~ /^$Echr1$/)&&
           ($pos2 <= $Epos2+$window)&&
            ($chr2 =~ /^$Echr2$/)&&
           ($pos3 <= $Epos3+$window)
          )
        {$nextline = <BEDPE> ;
		last if (!$nextline);
		$nextline=~chomp;
         ($chr1,$pos1,$pos2,$chr2,$pos3,$pos4)=split("\t",$nextline);
		 #print "NEXTLINE=$nextline";
         push (@chr1,$chr1);
         push (@pos1,$pos1);
         push (@pos2,$pos2);
         push (@chr2,$chr2);
         push (@pos3,$pos3);
         push (@pos4,$pos4);   
		  }
    ($Echr1,$Epos1,$Epos2,$Echr2,$Epos3,$Epos4)=($chr1[0],min(@pos1),max(@pos2),$chr2[-2],min(@pos3),$pos4[-2]);
    #print join("\t",$Echr1,$Epos1,$Epos2,$Echr2,$Epos3,$Epos4);
	if($pos1>$pos2){my $tmp=$pos1;$pos1=$pos2;$pos2=$tmp}
	if($pos1>$pos2){my $tmp=$pos3;$pos3=$pos4;$pos4=$tmp}
	print OUT1 join ("\t",$chr1,$pos1,$pos2)."\n";
	print OUT2 join ("\t",$chr2,$pos3,$pos4);
	($Echr1,$Epos1,$Epos2,$Echr2,$Epos3,$Epos4)=($chr1,$pos1,$pos2,$chr2,$pos3,$pos4);
	}
close BEDPE;
close OUT;
unlink ($tmpSpace2);

#####################################################################
#Now find out for each Unique BEDPE, how many Samples was the SV in?
#####################################################################
#FOR EACH VCF
#get NAME

my $tmpSpace5=join("",$tmpSpace,".5");
my $tmpSpace6=join("",$tmpSpace,".6");
my $tmpSpace7=join("",$tmpSpace,".7");
my $tmpSpace8=join("",$tmpSpace,".8");
my $tmpSpace9=join("",$tmpSpace,".9");

#Create a placeholder file
system("paste $tmpSpace3 $tmpSpace4| awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,0,\"NA\"}' > $tmpSpace7");
#Convert the VCF into a BED PE
for (my $i=0;$i<@VCF;$i++){
	open (OUT,">$tmpSpace5") or die "Cant write in this directory\n";
	open(VCF,$VCF[$i]) ;
	print STDERR "Starting on $VCF[$i]\n";
		while (<VCF>) {
			next if ($_=~/^#/);
			chomp;
			@line=split("\t",$_);
			$mate=$line[4];
			$mate=~s/[A-L]|[N-W]|[Z]|\[|\]//g;
			@mate=split(/:/,$mate);
			$end1a=$line[1]-$window;
			$end1b=$line[1]+$window;
			$end2a=$mate[1]-$window;
			$end2b=$mate[1]+$window;
			next if (($end1a<0)||($end2a<0));
			if (($line[0]=~/^chr$/)||($mate[0]=~/^chr$/)) {
				#print "$_\n";
				next;
			}
			print OUT "$line[0]\t$end1a\t$end1b\t$mate[0]\t$end2a\t$end2b\n";
			print OUT "$mate[0]\t$end2a\t$end2b\t$line[0]\t$end1a\t$end1b\n";
		}
	close VCF;
	close OUT;
	#for each row in $tmpSpace3, count the number of overlaps on both sides
	my $left=join("",$tmpSpace,".left");
	my $right=join("",$tmpSpace,".right");
	system("intersectBed -a $tmpSpace3 -b $tmpSpace5 -loj -c > $left");
	system("intersectBed -a $tmpSpace4 -b $tmpSpace5 -loj -c > $right");

	my $Lcount=`wc -l $left|cut -f1 -d" "`;
	my $Rcount=`wc -l $right|cut -f1 -d" "`;
	if ($Lcount != $Rcount){die "Need to check for errors in $left and $right\n\n"}
	system("paste $left $right > $tmpSpace5");
	system ("rm $left $right");
	open (IN,"$tmpSpace5") or die "Cant find $tmpSpace5\n";
	open (OUT,">$tmpSpace6") or die "Cant write in this directory\n";
	while(<IN>){
		$_=~chomp;
		@lines=split("\t",$_);
		if(($lines[3] > 0)&&($lines[6] > 0)){print OUT "1\t$VCF[$i]\n"}else{print OUT "0\t.\n"}
		}
	close IN;
	close OUT;

	system("paste $tmpSpace7 $tmpSpace6 > $tmpSpace8");
	#system("head $tmpSpace7 $tmpSpace8");
	 open (IN,"$tmpSpace8") or die "Cant find $tmpSpace8\n";
	 open (OUT,">$tmpSpace9") or die "Cant write in this directory\n";
	 my ($Samples,$NumSamples,$EVENT);
	 while(<IN>){
		 $_=~chomp;
		 @lines=split("\t",$_);

		 if ($lines[8] > 0){
			$Samples=$lines[7].";".$lines[9];
			$Samples=~s/^NA;//;
			$NumSamples=$lines[6]+$lines[8];
			}
			else{
			$Samples=$lines[7];
			$NumSamples=$lines[6];
			}
			print OUT join ("\t",@lines[0..5],$NumSamples,$Samples)."\n";
	 }
	 close IN;
	 close OUT;
	 print STDERR "completed with $VCF[$i]\n";
	 system("cp $tmpSpace9 $tmpSpace7");
}

system("cp $tmpSpace7 $outFile");
unlink ($tmpSpace9, $tmpSpace8, $tmpSpace7, $tmpSpace9,$tmpSpace3, $tmpSpace4, $tmpSpace5, $tmpSpace6);
print STDERR "Your results are in $outFile\n";


sub usage(){
    print "
###
### This script will merge multiple SoftSearch VCF files
###

Usage: Merge_SV.pl -v <vcf1> <vcf2> <vcfN> -w [500] -o <output file>
   
    Note: Must have bedtools installed and in your path\n\n\n";
}
