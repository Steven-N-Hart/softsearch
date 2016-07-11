#!/usr/bin/perl
#Author Steven Hart, PhD
#11-15-2012
#Convert and filter BAM files into merged bed 
#Output should be 
#chrA StartA EndA chrB StartB EndB Gene_id #supportingReads StrandA StrandB
#chr9 1000 5000 chr9 3000 3800 bedpe_example2 100 - +

use Cwd;
use File::Basename;
#Usage
sub usage(){
	print "Usage: perl Bam2Pair.pl -b <BAM> -o <outfile>\n
		-isize [10000]\t\tThe insert size to be considered discordant\n
		-winsize [10000]\tThe distance between mate pairs to be considered the same\n
		-min [1]\t\tThe minimum number of reads required to support an SV event\n
		-prefix need a random prefix so files with the same name don't get created\n\n"
		;
}
$bedtools=`which intersectBed`;
$samtools=`which samtools`;

if(!defined($bedtools)){die "BEDtools must be installed\n";}
if(!defined($samtools)){die "Samtools must be installed\n";}
use Getopt::Long;
#Declare variables
GetOptions(
	'b=s' => \$BAM_FILE,		#path to bam
	'out=s' => \$outfile,		#path to output
	'java:s' => \$java	,
        'chrom:s' => \$chrom      ,
	'isize=i' => \$isize,
	'winsize=i' => \$winsize,
        'prefix=s' => \$prefix,
	'min=i' => \$minSupport,
	'blacklist:s' => \$new_blacklist,
	'q=s' => \$qual,
	'v' => \$verbose
	);
#if(defined($picard_path)){$picard_path=$picard_path} else {die "Must specify a path to PICARD so that files can be sorted and indexed properly\n"};
if(!defined($BAM_FILE)){die "Must specify a BAM file!\n".usage();}
if(!defined($outfile)){die "Must specify an out filename!\n".usage();}
if(!defined($java)){$java=$java;}else{$java=`which java`}
if(!defined($qual)){$qual=20}
if($new_blacklist){$new_blacklist=" -L $new_blacklist"}


$Filter_BAM=$BAM_FILE;

@bam=split("/",$Filter_BAM);
$Filter_BAM=@bam[@bam-1];
$Filter_BAM=~s/.bam/$prefix.bam/;
$Filter_sam=$Filter_BAM;
$Filter_sam=~s/.bam/.sam/;




print "\nLooking for Discordant read pairs (and Unmated reads) without soft-clips\n";

#$command=join("","samtools view -h -q 20 -f 1 -F 1804 ",$BAM_FILE," ",$chrom," ",$new_blacklist," |  awk -F\'\\t\' \'{if (\$9 > ", $isize, " || \$9 < -",$isize," || \$9 == 0 || \$1 ~ /^@/) print \$0}' > ",$Filter_sam);


#Change command to allow reads where mate is unmapped & remove qual
$command=join("","samtools view -h -f 1 -F 1800 -q $qual ",$BAM_FILE," ",$chrom," ",$new_blacklist," |  awk -F\'\\t\' \'{if (\$9 > ", $isize, " || \$9 < -",$isize," || \$9 == 0 || \$1 ~ /^@/) print \$0}' > ",$Filter_sam);

print "$command\n" if ($verbose);
system($command);
$path = dirname(__FILE__);
$Filter_cluster=$Filter_sam;
$Filter_cluster=~s/.sam/.cluster/;
$command=join("",$path,"/ReadCluster.pl -i=$Filter_sam -o=$Filter_cluster -m $minSupport");
if($verbose){print "\n$command\n"};	

system($command);

##################################
#Now there are 2 SAM files of filtered reads
#.filter.cluster.inter.sam
#.filter.cluster.intra.sam
$result_pe=join("",$Filter_cluster,".out");
$command=join("","cat ",$Filter_cluster,".int\*|perl -ane 'next if(\@F[0]=~/^\@/);if(\@F[6]!~/=/){print join(\"\\t\",\$F[11],\@F[2],\@F[3],\@F[6],\@F[7],\"\\n\")}else{print join(\"\\t\",\$F[11],\@F[2],\@F[3],\@F[2],\@F[7],\"\\n\")}' >",$result_pe);
if($verbose){print $command."\n"};
system($command);
 #my ($sample, $chrstart, $start, $chrend, $end) 
$command=join("","cat ",$result_pe," | ",$path,"/cluster.pair.pl ",$winsize," |awk '(\$6 >",$minSupport,")' >> ", $outfile);
if($verbose){print $command."\n"};
system($command);
$filt1=join("",$Filter_cluster,".inter.sam");
$filt2=join("",$Filter_cluster,".intra.sam");


unlink($Filter_sam,$filt1,$filt2,$result_pe);

#########################################
#Now determine if left or righ clipping surrogate
print "\nBam2pair.pl Done\n";

