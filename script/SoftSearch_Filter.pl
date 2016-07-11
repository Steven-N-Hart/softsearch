#!/usr/bin/perl -s
open (FILE,"$ARGV[0]")||usage();#die "Not using the right Parameters!\n\n";
use Getopt::Long;
#Declare variables
my ($lsc,$minDist,$skip,$nSC,$nRP,$isize,$answer);
GetOptions(
	'dist:s' => \$minDist,		#minimum distance between events
	'lsc:i' => \$lsc,		#minimum somatic score
	'nsc:i' => \$nsc, 	#minimum depth of coverage in normal
	'nRP:i' => \$nRP,	#minimum number of times it can be seen in tumor
	'isize:i' => \$isize,	
	'sv:s' => \$sv,		#whether or not to skip small deletions
	'q:s' => \$answer,		#useful for plotting histograms
	'skip:s' => \$skip
	);
if(defined($lsc)){$lsc=$lsc} else {$lsc=0};
if(defined($nsc)){$nsc=$nsc} else {$nsc=0};
if(defined($nRP)){$nRP=$nRP} else {$nRP=0};
if(defined($minDist)){$minDist=$minDist} else {$minDist=0};
if(!$isize){$isize=0};
if(!$uRP){$uRP=0};

if($answer eq "yes"){$answer=$answer} else {$answer="no"};

if ($answer eq "yes"){
open(lsc,">lsc.out")||die;
open(nsc,">nsc.out")||die;
open(nRP,">nRP.out")||die;
}


#Remove hits if they are within $minDist
$chr="chr1";$pos=0;
while (<FILE>){
	if ($_=~/^#/){
		print; 
		next
	};
	if ($skip){next if $_=~/$skip/}
	@_=split(/\t/,$_);
	#Get ISIZE from INFO field
	my @info=split(/;/,$_[7]);
       	my $k = 0;
	my $v = 0; 
	my $infoHash;
	for (my $i=0;$i<=@info;$i++){
        	my @tmp=split(/=/,$info[$i]);
		$k=shift(@tmp);
		$v=shift(@tmp);
		$infoHash{$k}=$v;
	}

	#Get the value of TYPE to find out how many reads support the event
        my $counts = {CTX => 0, DEL => 0, INS => 0, INV => 0, TDUP => 0, NOV_INS => 0, lSC => 0, nSC => 0,uRP =>0,sDEL => 0,levD_local=>0,distl_levD => 0 };
	#Get Complete Hash
	#@_[8] is format
	#@_[9] is values
	my @format=split(/:/, $_[8]);
	my @sample=split(/:/,$_[9]);
	my %hash; 
	@hash{@format}=@sample;
	#Subset has to get proper type of variants
	my $max_val = 0;
	my $max_type = "NA";
	
	#Get TYPEOF HASH 
	my %type;
	%type = %hash ;
	delete $type{'lSC'};
        delete $type{'nSC'};
        delete $type{'uRP'};
        delete $type{'levD_local'};
        delete $type{'distl_levD'};

 	while (my ($key,$val)=each(%type)){
		if($val > $max_val){$max_val=$val;$max_type=$key}
		}


#######################################################################################################
        #Start applying filters
	
	#Remove hits if they are within $minDist
	$chrom=$_[0];$position=$_[1];

	#next if chroms are same and distance is less than X
	$difference=abs($pos-$position);
	if(($chrom eq $chr)&&($difference < $minDist)){
		$pos=$position;$chr=$chrom;;
		next}
	$pos=$position;$chr=$chrom;	
	$EVENT_SIZE=$infoHash{'ISIZE'};
	$EVENT_TYPE=$max_type;
	$EVENT_SUPPORT=$max_val;
	$length_of_softClips=$hash{'lSC'};
	$number_of_softclips=$hash{'nSC'};
        $number_of_unmated=$hash{'uRP'};
	
	########################################################################
	#Print if all fileds are ok
	next if($EVENT_SIZE < $isize);
        next if($EVENT_SUPPORT < $nRP);
        next if($length_of_softClips < $lsc);
        next if($number_of_softclips < $nsc);
        next if($number_of_unmated < $uRP);
	next if (($sv)&&($EVENT_TYPE=~/sDEL/));
	print;


	if ($answer eq "yes"){
	print lsc $length_of_softClips."\n";
	print nsc $number_of_softclips."\n";
	print nRP $EVENT_SUPPORT."\n";
	}
}


sub usage{
print "\nUsage: Soft_SearchFilter.pl <VCF>\n
	-dist	#minimum distance between events [0]
	-lsc	#minimum length soft-clip [0]
	-nsc	#minimum number of soft-clip [0]
	-nRP	#minimum number of discordant read pairs [0]
	-isize	#minimum size [0]
	-sv	#skip small deletions [no|yes]
	-skip	#pipe-delimited list of strings to skip (e.g. chrM|chY|chrGL)
	\n"
}

#R
# lsc<-read.table("lsc.out")
# nsc<-read.table("nsc.out")
# nRP<-read.table("nRP.out")
# par(mfrow=c(2,2))
# hist(lsc$V1,breaks=100)
# hist(nsc$V1,breaks=100)
# hist(nRP$V1,breaks=100)
