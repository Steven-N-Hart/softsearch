#!/usr/bin/perl -s
#Merge Softsearch results by chrom
if(!$ARGV[0]){die "Usage: <Sample.1.vcf>\n";}
my ($sample,$cmd);

#Get basename
$sample="$ARGV[0]";

$sample=~s/.[0-9(+)].out.vcf//;
#$sample=~s/.[0-9(+)].pe.vcf//;
#print "SAMPLE=$sample\n";
#exit;

my $outfile=$sample;
$outfile.="out.vcf";
if( -e $outfile ){unlink($outfile)}
$cmd="ls $sample\*vcf";
#print "$cmd\n";
my @samples=`$cmd`;
print "there are " .scalar(@samples)." samples\n";
#exit;
open (OUT,">$outfile");
my $i=1;
my $tmp=@samples[$i];
open(TMP,"$tmp");
while (<TMP>){
	print OUT if ($_=~/^#/);
}

open (OUT,">>$outfile");
my $chr;
for (my $i=0;$i<@samples;$i++){
	my $tmp=@samples[$i];
	open(TMP,"$tmp");
	while (<TMP>){
		unless (($_=~/^chrGL/)||($_=~/^#/)){print OUT $_;}
	}
	print "Done with $tmp";
        unlink($tmp);
	system("rm $tmp");
	close TMP;
}
close OUT;
