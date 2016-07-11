#!/usr/bin/perl

####
#### Usage: SoftSearch.pl [-lqrmsd] -b <BAM> -f <Genome.fa> -sam <samtools path> -bed <bedtools path>
#### Created 1-30-2012 by Steven Hart, PhD
#### hart.steven@mayo.edu
#### Required bedtools & samtools to be in path


use lib "/data2/bsi/reference/softsearch/lib" ;

use Getopt::Long;
use strict;
use warnings;
#use Data::Dumper;
use LevD;
use File::Basename;

my ($INPUT_BAM,$INPUT_FASTA,$OUTPUT_FILE,$minSoft,$minSoftReads,$dist_To_Soft,$bedtools,$samtools);
my ($minRP, $temp_output, $num_sd, $MapQ, $chrom, $unmated_pairs, $minBQ, $pair_only, $disable_RP_only);
my ($levD_local_threshold, $levD_distl_threshold,$pe_upper_limit,$high_qual,$sv_only,$blacklist,$genome_file,$verbose);

my $cmd = "";

#Declare variables
GetOptions(
	'b=s' => \$INPUT_BAM,
	'f=s' => \$INPUT_FASTA,
	'o:s' => \$OUTPUT_FILE,
	'm:i' => \$minRP,
	'l:i' => \$minSoft,
	'r:i' => \$minSoftReads,
	't:i' => \$temp_output,
	's:s' => \$num_sd,
	'd:i' => \$dist_To_Soft,
	'q:i' => \$MapQ,
	'c:s' => \$chrom,
	'u:s' => \$unmated_pairs,
	'x:s' => \$minBQ,
	'p' => \$pair_only,
	'g' => \$disable_RP_only,
	'j:s' => \$levD_local_threshold,
	'k:s' => \$levD_distl_threshold,
        'a:s' => \$pe_upper_limit,
        'e:s' => \$high_qual,
	'L' => \$sv_only,
	'v' => \$verbose, 
	'blacklist:s' => \$blacklist,
	'genome_file:s' => \$genome_file,
	"help|h|?"	=> \&usage);

unless($sv_only){$sv_only=""};
if(defined($INPUT_BAM)){$INPUT_BAM=$INPUT_BAM} else {print usage();die "Where is the BAM file?\n\n"}
if(defined($INPUT_FASTA)){$INPUT_FASTA=$INPUT_FASTA} else {print usage();die "Where is the fasta file?\n\n"}
my ($fn,$pathname) = fileparse($INPUT_BAM,".bam");
my $index=`ls $pathname/$fn*bai|head -1`;
#my $index =`ls \${INPUT_BAM%.bam}*bai`;
#print "INDEX=$index\n";
if(!$index){die "\n\nERROR: you need index your BAM file\n\n"}

### get current time
print "Start Time : " . &spGetCurDateTime() . "\n";
my $now = time;

#if(defined($OUTPUT_FILE)){$OUTPUT_FILE=$OUTPUT_FILE} else {$OUTPUT_FILE="output.vcf"; print "\nNo outfile specified.  Using output.vcf as default\n\n"}
if(defined($minSoft)){$minSoft=$minSoft} else {$minSoft=5}
if(defined($minRP)){$minRP=$minRP} else {$minRP=5}
if(defined($minSoftReads)){$minSoftReads=$minSoftReads} else {$minSoftReads=5}
#if(defined($dist_To_Soft)){$dist_To_Soft=$dist_To_Soft} else {$dist_To_Soft=300} #Changed to max normal insert size
if(defined($num_sd)){$num_sd=$num_sd} else {$num_sd=6}
if(defined($MapQ)){$MapQ=$MapQ} else {$MapQ=20}

unless (defined $pe_upper_limit) { $pe_upper_limit = 10000; }
unless (defined $levD_local_threshold) { $levD_local_threshold = 0.05; }
unless (defined $levD_distl_threshold) { $levD_distl_threshold = 0.05; }
#Get sample name if available
my $SAMPLE_NAME="";
my $OUTNAME ="";
$SAMPLE_NAME=`samtools view -f2 -H $INPUT_BAM|awk '{if(\$1~/^\@RG/){sub("ID:","",\$2);name=\$2;print name}}'|head -1`;
$SAMPLE_NAME=~s/\n//g;
if (!$OUTPUT_FILE){
	if($SAMPLE_NAME ne ""){$OUTNAME=$SAMPLE_NAME.".vcf"}
	else {$OUTNAME="output.vcf"}
}
else{$OUTNAME=$OUTPUT_FILE}

print "Writing results to $OUTNAME\n";


##Make sure if submitting on SGE, to prepned the "chr".  Not all referecne FAST files require "chr", so we shouldn't force the issue.
if(!defined($chrom)){$chrom=""}
if(!defined($unmated_pairs)){$unmated_pairs=0}

my $badQualValue=chr($MapQ);
if(defined($minBQ)){ $badQualValue=chr($minBQ); }

if($badQualValue  eq "#"){$badQualValue="\#"}

# adding and cheking for samtools and bedtools in the PATh
## check for bedtools and samtools in the path
$bedtools=`which intersectBed` ;
if(!defined($bedtools)){die "\nError:\n\tno bedtools. Please install bedtools and add to the path\n";}
#$samtools=`samtools 2>&1`;
$samtools=`which samtools`;
if($samtools !~ /(samtools)/i){die "\nError:\n\tno samtools. Please install samtools and add to the path\n";}

print "Usage = SoftSearch.pl -l $minSoft -q $MapQ -r $minSoftReads -d $dist_To_Soft -m $minRP -s $num_sd -c $chrom -b $INPUT_BAM -f $INPUT_FASTA -o $OUTNAME \n\n";
sub usage {
	print "\nusage: SoftSearch.pl [-cqlrmsd] -b <BAM> -f <Genome.fa> \n";
	print "\t-q\t\tMinimum mapping quality [20]\n";
	print "\t-l\t\tMinimum length of soft-clipped segment [5]\n";
	print "\t-r\t\tMinimum depth of soft-clipped reads at position [5]\n";
	print "\t-m\t\tMinimum number of discordant read pairs [5]\n";
	print "\t-s\t\tNumber of sd away from mean to be considered discordant [6]\n";
	print "\t-u\t\tNumber of unmated pairs [0]\n";
	print "\t-d\t\tMax distance between soft-clipped segments and discordant read pairs [Maximum normal insert]\n";
	print "\t-o\t\tOutput file name [output.vcf]\n";
	print "\t-t\t\tPrint temp files for debugging [no|yes]\n";
	print "\t-c\t\tuse only this chrom or chr:pos1-pos2\n";
	print "\t-p\t\tuse paired-end mode only. In other words, don't try to find soft-clipping events!\n";
	print "\t-g\t\tEnable paired-only seach. This will look for discordant read pairs even without soft clips.\n";
        print "\t-a\t\tset the minimum distance for a discordant read pair without soft-clipping info [10000]\n";
        print "\t-L\t\tFlag to print out even small deletions (low quality)\n";
        print "\t-e\t\tdisable strict quality filtering of base qualities in soft-clipped reads [no]\n";
        print "\t-blacklist\tareas of the genome to skip calling.  Requires -genome_file\n";
        print "\t-genome_file\ttab seperated value of chromosome name and length.  Only used with -blacklist option\n\n";

	exit 1;
	}


#############################################################
# create temporary variable name
#############################################################
srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
our $random_name=join "", map { ("a".."z")[rand 26] } 1..8;

#############################################################
## create green list
##############################################################
#
my $new_blacklist="";
if($blacklist){
        if(!$genome_file){die "if using a blacklist, you must also specify the location of a genome_file
        The format of the genome_file should be
                chrom   size
                chr1    249250621
                chr2    243199373
                ...

        If using hg19, you can ge the genome file by
                mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from hg19.chromInfo\"  > hg19.genome";}
        
	$cmd=join("","complementBed -i $blacklist -g $genome_file >",$random_name,".bed") ;
	system ($cmd);
	$new_blacklist=join(""," -L ",$random_name,".bed ");
	}

if($verbose){print "CMD=$cmd\nBlacklist is $new_blacklist\n";}





#############################################################
# Calcualte insert size distribution of properly mated reads
#############################################################

#Change for compatability with other operating systems
#my $metrics=`samtools view -q $MapQ -f2 $INPUT_BAM $chrom|cut -f9|head -10000|awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}'`;

my $metrics=`samtools view -q $MapQ -f2 $INPUT_BAM $chrom|cut -f9|head -10000|awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}'`;
#my ($mean,$stdev)=split(/ /,$metrics);
my ($mean,$stdev)=split(/\s/,$metrics);
$stdev=~s/\n//;
my $upper_limit=int($mean+($num_sd*$stdev));
my $lower_limit=int($mean-($num_sd*$stdev));
die if (!$mean);
print qq{The mean insert size is $mean +/- $stdev (sd)
The upper limit = $upper_limit
The lower limit = $lower_limit\n
};
if($lower_limit<0){
	print "Warning!! Given this insert size distribution, we can not call small indels.  No other data will be affected\n";
	$lower_limit=1;
}
my $tmp_name=join ("",$random_name,".tmp.bam");
my $random_file_sc = "";
my $command = "";
if(defined($dist_To_Soft)){$dist_To_Soft=$dist_To_Soft} else {$dist_To_Soft=$upper_limit}

#############################################################
# Make sam file that has soft clipped reads
#############################################################
#give file a name
if(!defined($pair_only)){
	$random_file_sc=join ("",$random_name,".sc.sam");
	$command=join ("","samtools view -q $MapQ -F 1024 $INPUT_BAM $chrom $new_blacklist| awk '{OFS=\"\\t\"}{c=0;if(\$6~/S/){++c};if(c == 1){print}}' | perl -ane '\$TR=(\@F[10]=~tr/\#//);if(\$TR<2){print}' > ", $random_file_sc);

	print "Making SAM file of soft-clipped reads\n";
if($verbose){	print "$command\n";}
	system("$command");

	#############################################################
	# Find areas that have deep enough soft-clip coverage
	print "Identifying soft-clipped regions that are at least $minSoft bp long \n";
	open (FILE,"$random_file_sc")||die "Can't open soft-clipped sam file $random_file_sc\n";

	my $tmpfile=join("",$random_file_sc,".sc.passfilter");
	open (OUT,">$tmpfile")||die "Can't write files here!\n";

	while(<FILE>){
		@_ = split(/\t/, $_);
		#### parse CIGAR string and create a hash of array of each operation
		my @CIGAR = split(/([0-9]+[SMIDNHXP])/, $_[5]);
		my $hash;
		map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;

		#for ($i=0; $i<=$#softclip_pos; $i++)	{
		foreach my $softclip (@{$hash->{S}}) {
			#if	($CIGAR[$softclip_pos[$i]] > $minSoft){
			if	($softclip > $minSoft){
				###############Make sure base qualities don't have more than 2 bad marks
				my $qual=$_[10];
				my $TR=($qual=~tr/$badQualValue//);
				if($badQualValue eq "#"){ $TR=($qual=~tr/\#//); }
				#Skip the soft clip if there is more than 2 bad qual values
				#next if($TR > 2);
#				if (!$high_qual){next if($TR > 2);}
				print OUT;
				last;
			}
		}
	}
	close FILE;
	close OUT;

	$command=join(" ","mv",$tmpfile,$random_file_sc);
if($verbose){	print "$command\n";}
	system("$command");
}

#########################################################
#Stack up SoftClips
#########################################################
my $random_file=join("",$random_name,".sc.direction.bed");
if(!defined($pair_only)){
        open (FILE,"$random_file_sc")|| die "Can't open sam file\n";
        #$random_file=join("",$random_name,".sc.direction");

        print "Calling sides of soft-clips\n";
        #\nTMPOUT=$random_file\tINPUT=$random_file_sc\n\n";
        open (TMPOUT,">$random_file")|| die "Can't create tmp file\n";

        while (<FILE>){
                @_ = split(/\t/, $_);
                #### parse CIGAR string and create a hash of array of each operation
                my @CIGAR = split(/([0-9]+[SMIDNHXP])/, $_[5]);
                my $hash;
                map { push(@{$hash->{$2}}, $1) if (/(\d+)([SMIDNHXP])/) } @CIGAR;

                #### next if softclips on each end
                next if ($_[5] =~ /^[0-9]+S.*S$/);

                #### next softclip occurs in the middle
                #next if ($_[5] =~ /^[0-9]+[^S][0-9].*S.+$/);
                #Thanks to Michael Ta for fixing this, now softclips > 100 bp will work
		next if ($_[5] =~ /^[0-9]+[A-Z][0-9].*S.+$/);
                my $softclip = $hash->{S}[0];

                my $end1 = 0;
                my $end2 = 0;
                my $softBases = "";
		my $right_corrected="";my $left_corrected="";

        if ($softclip > $minSoft) {
		
                        ####If the soft clip occurs at end of read and its on the minus strand, then it's a right clip
                        if ($_[5] =~ /^.*S$/) {
                                $end1=$_[3]+length($_[9])-$softclip-1;
                                $end2=$end1+1;
next if ($end1<0);
                                #RIGHT clip on Minus
                                $softBases=substr($_[9], length($_[9])-$softclip, length($_[9]));
                                #Right clips don't always get clipped correctly, so fix that
                                # Check to see if sc base matches ref
                                $right_corrected=baseCheck($_[2],$end2,"right",$softBases);
                               print TMPOUT "$right_corrected\n"

                        } else {
                                #### Begins with S (left clip)
                                $end1=$_[3]-$softclip;
next if ($end1<0);

                                $softBases=substr($_[9], 0,$softclip);#print "TMP=$softBases\n";
        			$left_corrected=baseCheck($_[2],$end1,"left",$softBases);
if(!$left_corrected){print "baseCheck($_[2],$end1,left,$softBases)\n";next}
                               print TMPOUT "$left_corrected\n";
#print "\nSEQ=$_[9]\t\n";

                        }
        }
  }
close FILE;
close TMPOUT;
}
sub baseCheck{
        my ($chrom,$pos,$direction,$softBases)=@_;
        #skip if position is less than 0, which is caused by MT DNA
        return if ($pos<0);
        my $exit="";

        while(!$exit){
        if($direction=~/right/){
                        my $refBase=getSeq($chrom,$pos,$INPUT_FASTA);
                        my $softBase=substr($softBases,0,1);
                        if ($softBase !~ /$refBase/){
                                my $value=join("\t",$chrom,$pos,$pos+1,join("|",$softBases,$direction));
                                $exit="STOP";
                                return $value;
                        }
                        else{
                                $pos=$pos+1;
                                $softBases=substr($softBases, 1,length($softBases));
                        }
         }
        else{
                        my $refBase=getSeq($chrom,$pos+1,$INPUT_FASTA);
                        my $softBase=substr($softBases,-1,1);
                        if ($softBase !~ /$refBase/){
                                $pos=$pos-1+length($softBases);
                                my $value=join("\t",$chrom,$pos-1,$pos,join("|",$softBases,$direction));
                                $exit="STOP";
                                return $value;
                        }
                        else{
                                $pos=$pos-1;
                                $softBases=substr($softBases, 0, -1);
                                #print "Trying again $softBases\n";
                       }

        }

}
}
#Remove SAM files to conserve space
unlink($random_file_sc);


my $random_file_disc="$INPUT_BAM";
###
#
######################################################
# Transform Read pair groups into softclip equivalents
######################################################
#
#
#
my $v="";
#if($disable_RP_only){
print "Running Bam2pair.pl\n";
print "Looking for discordant read pairs without requiring soft-clipping information\n";
	use FindBin qw($Bin);
	my $path=$Bin;
#	print"\n\nPATH=$path\n\n";
if($verbose){$v="-v"}
	my $tmp_out=join("",$random_name,"RP.out");
##Provided by user
	if ($new_blacklist ne "") {
		$command=join("","perl ",$path,"/Bam2pair.pl -b $random_file_disc  -o $tmp_out -isize $pe_upper_limit -winsize $dist_To_Soft -min $minRP -chrom $chrom -prefix $random_name -q $MapQ -blacklist $random_name.bed $v");
	} else {
		$command=join("","perl ",$path,"/Bam2pair.pl -b $random_file_disc  -o $tmp_out -isize $pe_upper_limit -winsize $dist_To_Soft -min $minRP -chrom $chrom -prefix $random_name -q $MapQ $v")
	}

if($verbose){	print "$command\n"};
	system("$command");
	$command=join("","perl -ane '\$end1=\@F[1];\$end2=\@F[3];print join(\"\\t\",\@F[0..1],\$end1,\"unknown|left\");print \"\\n\";print join(\"\\t\",\@F[2..3],\$end2,\"unknown|left\");print \"\\n\"' ", $tmp_out," >> ",$random_file);
if($verbose){print "$command\n"};
	system($command);
	unlink($tmp_out);
#}
#


######################################################
unlink("$random_file","$tmp_name","$random_file","$index","$random_name","$new_blacklist") if (-z $random_file || ! -e $random_file ) ;
if (-z $random_file || ! -e $random_file){
	print "Softclipped file is empty($random_file).\nNo soft clipping found using desired paramters\n\n";
	open (OUT,">$OUTNAME")||die "Can't write files here!\n";
        &print_header();
        close OUT;
        exit;
	}


#############################################################
#  Make sure there are enough soft-clippped supporting reads
#############################################################
my $outfile=join("",$random_file,".sc.merge.bed");
#sortbed -i .sc.direction | mergeBed -nms -d 25 -i stdin > .sc.merge.bed
$command=join(" ","sortBed -i",$random_file," | mergeBed  -nms -i stdin","|egrep \";|,\"","|awk '{OFS=\"\t\"}(NF==4)'",">",$outfile);

print "$command\n" if ($verbose);
system("$command");

if (-z $outfile || ! -e $outfile){
	unlink("$tmp_name","$random_file","$outfile","$index","$random_name","$new_blacklist"); 
	print "mergeBed file is empty.\nNo strucutral variants found\n\n" ;
        open (OUT,">$OUTNAME")||die "Can't write files here!\n";
        &print_header();
        close OUT;
        exit;
}

print "completed mergeBed\n";

###############################################################
# If left and right are on the same line, make into 2 lines
###############################################################
open (INFILE,$outfile)||die "couldn't open temp file : $. \n\n";
my $tmp2=join("",$random_name,".sc.fixed.merge.bed");
#print "INFILE=$outfile\tOUTFILE=$tmp2\n\n";
#INPUT FORMAT=chr9\t131467\t131473\tATGCTTATTAAAA|left;TTATTAAAAGCATA|left
open (OUTFILE,">$tmp2")||die "couldn't create temp file : $. \n\n";
while(<INFILE>){
	chomp $_;
	my $l = $_;

	my @a = split(/\t/, $l);
	my $info = $a[3];
	my @info_arr = split(/\;/, $info);
	my @left_arr=();
	my @right_arr=();
	@left_arr = grep(/left/, @info_arr);
	@right_arr = grep(/right/, @info_arr);

	#New
	my $left = join(";", @left_arr);
	my $right = join(";", @right_arr);
	$info = join(";", @info_arr);

	if((@left_arr) && (@right_arr)){
		print OUTFILE "$a[0]\t$a[1]\t$a[2]\t$left\n$a[0]\t$a[1]\t$a[2]\t$right\n";
	} else{
		my $all=join("\t",@a[0..2],$info);
		print OUTFILE "$all\n";
	}
}

# make sure output file name is $outfile
$command=join(" ","sed -e '/ /s//\t/g'", $tmp2,"|awk 'BEGIN{OFS=\"\\t\"}(NF==4)'", "|perl -pne 's/ /\t/g'>",$outfile);
system("$command");
if($verbose){print "$command\n"};
unlink("$tmp_name","$random_file","$tmp2","$outfile","$index","random_name","$new_blacklist") if (-z $outfile || ! -e $outfile) ;
 if (-z $outfile || ! -e $outfile){
	print "Fixed mergeBed file is empty($outfile).\nNo strucutral variants found\n\n";
        open (OUT,">$OUTNAME")||die "Can't write files here!\n";
        &print_header();
        close OUT;
        exit;
}

print "completed fixing mergeBed\n\n";

###############################################################
# Seperate directions of soft clips
###############################################################
my $left_sc = join("", "left", $tmp2);
my $right_sc = join("", "right", $tmp2);
use FindBin qw($Bin);
#my $path=$Bin;

$command=join("","grep left ", $tmp2, " |sed -e '/left /s//left\;/g' | sed -e '/ /s//\t/g'|perl ".$path."/direction_filter.pl - >",$left_sc);
system("$command");
#print "$command\n";
$command=join("","grep right ", $tmp2, " |sed -e '/right /s//right\;/g' | sed -e '/ /s//\t/g'|perl ".$path."/direction_filter.pl - >",$right_sc);
#$command=join(" ","grep right ", $tmp2, " |sed -e '/right /s//right\;/g' | sed -e '/ /s//\t/g' >",$right_sc);
system("$command");
#print "$command\n";
#die "CHECK $right_sc\n";

###############################################################
# Count the number and identify directions of soft clips
###############################################################
print "Count the number and identify directions of soft clips\n";
#print "looking in $outfile\n";
$outfile=join("",$random_name,".sc.fixed.merge.bed");

open (INFILE,$outfile)||die "couldn't open temp file\n\n";
my $tmp3 = join("", $random_file, "predSV");
open (OUTFILE, ">$tmp3")||die "couldn't create temp file\n\n";
while(<INFILE>){
chomp;
	@_=split(/\t/,$_);
	my $count=tr/\;//;$count+=tr/\,//;
	$count=$count+1;
	my $left=0;
	my $right=0;

	while ($_ =~ /left/g) { $left++ } # count number of right clips
	while ($_ =~ /right/g) { $right++ } # count number of left clips

	###############################################################
	if ($count >= $minSoftReads){
		####get longets soft-clipped read
		my @clips=split(/\;|,|\|/,$_[3]);

		my ($max, $temp, $temp2, $temp3, $dir, $maxSclip) = (0) x 6;
		for (my $i=0; $i<$count; $i++) {
			my $plus1=$i+1;
			$temp=length($clips[$i]);
			$temp2=$clips[$plus1];
			$temp3=$clips[$i];

			if ($temp > $max){
				$maxSclip=$temp3;
				$max =$temp;
				$dir=$temp2;
			} else {
				$max=$max;
				$dir=$dir;
				$maxSclip=$maxSclip;
			}
			$i++;
		}
		my $order2 = join("|", $left, $right);
        #print join ("\t",@_[0..2],$maxSclip,$max,$dir,$count,$order2) . "\n";
		print OUTFILE join ("\t",@_[0..2],$maxSclip,$max,$dir,$count,$order2) . "\n";
	} elsif($_=~/unknown/){
	print OUTFILE join ("\t",@_[0..2],"NA","NA","left","NA","NA|NA") . "\n";
        print OUTFILE join ("\t",@_[0..2],"NA","NA","right","NA","NA|NA") . "\n";
	}
	####Format is Chrom,start, end,longest Soft-clip,length of longest Soft-clip, direction of longest soft-clip,#supporting softclips,#right Sclips|#left Sclips
}
close INFILE;
close OUTFILE;

unlink("$tmp2","$tmp_name","$random_file","$tmp3","$outfile","$index","$random_name","$right_sc","$left_sc","$new_blacklist") if (-z $tmp3 || !-e $tmp3) ;

 if (-z $tmp3 || !-e $tmp3){
	print "No structural variants found while Counting the number and identify directions of soft clips.\n" ;

	open (OUT,">$OUTNAME")||die "Can't write files here!\n";
	&print_header();
	close OUT;
	exit;

}

print "Done counting Softclipped reads\n";
###############################################################
#### Print header information
###############################################################
open (OUT,">$OUTNAME")||die "Can't write files here!\n";
&print_header();
close OUT;

###############################################################
###############################################################
#### DO the bulk of the work
###############################################################
use List::Util qw(min max);
open (FILE,"$tmp3")|| die "Can't open file\n";
open (OUT,">>$OUTNAME")|| die "Can't open file\n";

#print "\nusing $tmp3 and writing to $OUTPUT_FILE \n";
while (<FILE>){
	#If left clip {+- or -- or -+ }{+- are uninformative b/c they go upstream}
	#If right clip {++ or -- or +-}
	chomp $_;
	my $line = $_;

	my @info = split(/\t/, $_);

#INFO = chr5 49559445 49559446 CCTACTATAGGCATGAAAGTGCTGCAAATATCCGTTT 37 right 5 0|5
#Count number of alt chroms
#Get bases
my $Start=$info[1]-1000 ;
my $End=$info[2]+1000 ;
my $res=`samtools view $INPUT_BAM $info[0]:$Start-$End|cut -f7|sort|uniq -c|wc -l|awk '{print \$1}'`;
next if ($res > 4);
	if($info[5] eq "left") {
		bulk_work("left", $line, $random_file_disc);

	} elsif ($info[5] eq "right") {
		bulk_work("right", $line, $random_file_disc);
	}
#if($. ==6){print "THIS IS LINE 6\n$_\n";die}
print "Completed line $.\n" if ($verbose);
}
close FILE;
close OUT;

###############################################################################
###############################################################################
#### Delete temp files
my $meregedBed=join("",$random_name,".sc.direction.bed.sc.merge.bed");

if(defined($temp_output)){$temp_output=$temp_output} else {$temp_output="no"}

if ($temp_output eq "no"){
	unlink("$tmp_name","$random_file","$tmp2",,"$tmp3","$outfile","$index","$random_name","$right_sc","$left_sc","$meregedBed","$random_name.bed");
}
####Sort VCF
my $tmp=join(".",$random_name,"tmp");
#Get header
$cmd="grep \"#\" $OUTNAME > $tmp";
system($cmd);
#sort results
$cmd="grep -v \"#\" $OUTNAME|perl -pne 's/chr//'|sort -k1,1n -k2,2n|perl -ne 'print \"chr\".\$_' >>$tmp";
system($cmd);
$cmd="mv $tmp $OUTNAME";
system($cmd);
#remove entries next to each other




#############################################################
##May not need this anymore since filtering on left and right
#############################################################
#my $tmpout=$OUTNAME;
#$tmpout.=".tmp";
#use FindBin qw($Bin);
##my $path=$Bin;
#$command="perl ".$path."/Extract_nSC.pl $OUTNAME -q nSC > $tmpout";
##print "Command=$command\n";
#system($command);
#$command="perl ".$path."/reduce_redundancy.pl $tmpout $upper_limit |cut -f1-10 > $OUTNAME";
##print "$command\n";
#system($command);
#system("rm $tmpout");
########################################################




print "Analysis Completed\n\nYou did it!!!\n";
print "Finish Time : " . &spGetCurDateTime() . "\n";
$now = time - $now;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60),
int($now % 60));

exit;

###############################################################################
sub rev_comp {
  my $dna = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}


###############################################################################
#### to get reference base
sub getSeq{
	my ($chr,$pos,$fasta)=@_;
	#don't require chr
	#if($chr !~ /^chr/){die "$chr is not correct\n";}
#	die "$pos is not a number\n" if ($pos <0);
my @result=();
        if ($pos <0){print "$pos is not a valid position (likely caused by circular MT chromosome)\n";return;}

	@result = `samtools faidx $fasta $chr:$pos-$pos`;
	if($result[1]){chomp($result[1]);
	return uc($result[1]);
	}
	return("NA");
	#### after return will not be printed
	####print "RESULTS=@result\n";
}

sub getBases{
        my ($chr,$pos1,$pos2,$fasta)=@_;
        #don't require chr
        #if($chr !~ /^chr/){die "$chr is not correct\n";}
my @result=();
        if ($pos1 <0){print "$pos1 is not a valid position (likely caused by circular MT chromosome)\n";return;};

        @result = `samtools faidx $fasta $chr:$pos1-$pos2`;
	if(!$result[1]){$result[1]="NA"};
        chomp($result[1]);
        return uc($result[1]);

        #### after return will not be printed
        ####print "RESULTS=@result\n";
}
###############################################################################
#### to get time
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}


###############################################################################
#### print header
sub print_header {
	my $date=&spGetCurDateTime();
	my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=SoftSearch.pl
##reference=$INPUT_FASTA
##Usage= SoftSearch.pl -l $minSoft -q $MapQ -r $minSoftReads -d $dist_To_Soft -m $minRP -u $unmated_pairs -s $num_sd -b $INPUT_BAM -f $INPUT_FASTA -o $OUTNAME
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=ISIZE,Number=.,Type=String,Description="Size of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##FORMAT=<ID=lSC,Number=1,Type=Integer,Description="Length of the longest soft clips supporting the BND">
##FORMAT=<ID=nSC,Number=1,Type=Integer,Description="Number of supporting soft-clips\">
##FORMAT=<ID=uRP,Number=1,Type=Integer,Description="Number of unmated read pairs nearby Soft-Clips">
##FORMAT=<ID=levD_local,Number=1,Type=Float,Description="Levenstein distance between soft-clipped bases and the area around the original soft-clipped site">
##FORMAT=<ID=distl_levD,Number=1,Type=Float,Description="Levenstein distance between the soft-clipped bases and mate location">
##FORMAT=<ID=CTX,Number=1,Type=Integer,Description="Number of chromosomal translocations">
##FORMAT=<ID=DEL,Number=1,Type=Integer,Description="Number of reads supporting Large Deletions">
##FORMAT=<ID=INS,Number=1,Type=Integer,Description="Number of reads supporting Large insertions">
##FORMAT=<ID=NOV_INS,Number=1,Type=Integer,Description="Number of reads supporting novel sequence insertion">
##FORMAT=<ID=TDUP,Number=1,Type=Integer,Description="Number of reads supporting a tandem duplication">
##FORMAT=<ID=INV,Number=1,Type=Integer,Description="Number of reads supporting inversions">
##FORMAT=<ID=sDEL,Number=1,Type=Integer,Description="Number of reads supporting novel sequence insertion">
##INFO=<ID=NO_MATE_SC,Number=1,Type=Flag,Description="When there is no softclipping of the mate read location, an appromiate position is used">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Dummy value for maintaining VCF-Spec">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$SAMPLE_NAME\n};

	print OUT $header;
}


###############################################################################
sub bulk_work {
print "#####################################@_\n" if ($verbose);
	my ($side, $line, $file) = @_;
	my $local_levD = 0;
	my $distl_levD = 0;

	#my @info = split(/\t/, $line);
	my @plus_Reads = split(/\t/, $line);
	$plus_Reads[7] =~ s/\n//g;

	#### softclip length and softclip size.
	my $lSC = $plus_Reads[4];
	my $nSC = $plus_Reads[6];


	#Get all types of compatible reads
	#Get improperly paired reads (@ max distance)

	#### default value for left SIDE.
	#If left-clip, then look downstream for match of softclipped reads to define a deletion, but look for DRPs upstream
	my $sv_type = "SVTYPE=BND";
	my $start_local=0; my $end_local=0;my $target_local="";my $target_drp="";my $start_drp="";my $end_drp="";
	if ($side =~ /left/) {
		$start_local = $plus_Reads[1]-$dist_To_Soft;
		$end_local = $plus_Reads[2];
                $start_drp = $plus_Reads[1];
                $end_drp = $plus_Reads[1]+$dist_To_Soft;
	
	}
	else{                
                $start_local = $plus_Reads[1];
                $end_local = $plus_Reads[1]+$dist_To_Soft;
                $start_drp = $plus_Reads[1]-$dist_To_Soft;
                $end_drp = $plus_Reads[1];
        }
	
	$target_local=join("", $plus_Reads[0], ":", $start_local, "-", $end_local);
	$target_drp=join("", $plus_Reads[0], ":", $start_drp, "-", $end_drp);
	my $num_unmapped_pairs="";
	if ($side =~ /right/) {
		$num_unmapped_pairs=`samtools view $new_blacklist -q $MapQ -f8 -F 1536 -c $INPUT_BAM $target_drp`;
	} else {
        $num_unmapped_pairs=`samtools view $new_blacklist -q $MapQ -f24 -F 1536 -c $INPUT_BAM $target_drp`;
	}
if($verbose){print "samtools view $new_blacklist -q $MapQ -f24 -F 1536 -c $INPUT_BAM $target_drp\n";}

	$num_unmapped_pairs=~s/\n//;
if($verbose){print "NUM UNMAPPED PAIRS= $num_unmapped_pairs\n";}
	my $REF1_base = "";
	my $REF2_base = "";
	my $INFO_1 = "";
	my $INFO_2 = "";
	my $ALT_1 = "";
	my $ALT_2 = "";
	my $isize = 0;
	my $QUAL = "";
	my $FORMAT = "GT:";

	#### get 8 bit rand id
	my $BND1_name = join "", map { ("a".."z")[rand 26] } 1..8;
	my $BND2_name = join "", map { ("a".."z")[rand 26] } 1..8;
	$BND1_name=join "_","BND",$BND1_name;
	$BND2_name=join "_","BND",$BND2_name;

	my $counts = {CTX => 0, DEL => 0, INS => 0, INV => 0, TDUP => 0, NOV_INS => 0 };
	my $event_mate_info = {CTX => "", DEL => "", INS => "", INV => "", TDUP => "", NOV_INS => "" };

	#### get mate pair info and counts per event
	foreach my $e (sort keys %{$counts}) {
		my $h = get_counts_n_info($e, $side, $MapQ, $file, $dist_To_Soft, $target_drp, $upper_limit, $lower_limit);

		$counts->{$e} = $h->{count};
		$event_mate_info->{$e} = $h->{info};
	}
#print Dumper($counts);

	my $max = 0;
	my $type = "UNKNOWN";
	my $nRP = 0;
	my $mate_info = "NA\tNA\tNA\tNA";
	my $summary = "GT:";

	#### find max count of events and set type, nRP and info to corresponding
	#### max count event.
	#### also create a summary string of all counts to be added to VCF file.
	foreach my $e (sort keys %{$counts}){
#		if ($counts->{$e} >=i $max){
		if ($counts->{$e} > $max){		
			$type = $e .",". $counts->{$e};
			$nRP = $counts->{$e};

			$max = $counts->{$e};

			if (length($event_mate_info->{$e})) {
				$mate_info = $event_mate_info->{$e};
			}
		}

		$summary .= $e .",". $counts->{$e} .":";
	}
#	print "done with Summary\n";
	#### remove last colon ":" from
	$summary =~ s/:$//;
 if (($minRP > $max)&&(!$disable_RP_only )){if ($verbose){print "FAILED BECAUSE ($minRP > $max)&&(!$disable_RP_only )"};return};

	#### Run Levenstein distance on softclip in target region to find out if its a small deletion/insetion
	#### passing 1: clip_seq, 2: chr, 3: start, 4: end, 5: ref file.
	my $levD = new LevD;
########################################################
########################################################
########################################################

	#### redefine start and end location for LevD calc.
#	$start = $plus_Reads[1]-$dist_To_Soft;
#	$end = $plus_Reads[2];
	my $num_bases_to_loc=0;
	my $new_start=0;
	my $new_end=0;
	my $del_seq="";
        my $start = $start_local;
        my $end = $end_local;
	if ($lSC=~/NA/){$lSC=0}

	if ($side =~ /right/) {
	        $levD->search($plus_Reads[3], $plus_Reads[0], $start, $end, $INPUT_FASTA);
		$local_levD = sprintf("%.2f", $levD->{relative_edit_dist});
	        $num_bases_to_loc=$levD->{index};
		$new_start = $plus_Reads[2];
                if ($plus_Reads[2]=~/^[0-9]/){$new_end=$plus_Reads[2]+$lSC};
	}
	else{
		$levD->search($plus_Reads[3], $plus_Reads[0], $start, $end, $INPUT_FASTA);
		$local_levD = sprintf("%.2f", $levD->{relative_edit_dist});
		$num_bases_to_loc=$levD->{index};
		if ($plus_Reads[2]=~/^[0-9]/){$new_start=$plus_Reads[2]-$lSC};
                $new_end = $plus_Reads[2];
	}
	if((!$new_start)||(!$new_end)||($new_start<0)){print "FAILED AT ((!$new_start)||(!$new_end)||($new_start<0))\n";return};
	
	$del_seq=getBases($plus_Reads[0], $new_start,$new_end,$INPUT_FASTA);
##############################################################################
#	#If there is a match, where is the start position of the match?
#
##############################################################################


	#if $plus_Reads[3] eq "NA", then it was found without soft-clipped reads
	if($plus_Reads[3] !~  /NA/){
			if (($local_levD < $levD_local_threshold)) {
				return if (!$sv_only);
				#### add value to summary to be written to vcf file.
				$summary = "GT:sDel," . $plus_Reads[6];
				$type = "sDEL";
				###########################################################################
				##### Printing output

				#########################################
				##### Get DNA info
				#########################################
				#$REF1_base = getSeq($plus_Reads[0], $plus_Reads[1], $INPUT_FASTA);
				$REF1_base = substr($del_seq, 0, 1);

				#### this is alt ref. for softclip its $plus_Reads[3]
				$REF2_base = $del_seq;
				$QUAL = 1/($local_levD + 0.001);
				$QUAL = sprintf("%.2f",$QUAL);
				$isize = length($del_seq);

				#### svtype = none for sDEL
				#### isize = length($info[3]);
				#### nRP = NA
				#### mate_id = NA
				#### CTX,:DEL,:....sDEL,##
				$INFO_1=join(";", "SVTYPE=NA", "EVENT=$type", "ISIZE=$isize");

				#Add Sample infomration
				my $FORMAT="GT:sDEL";
				$FORMAT .= ":lSC:nSC:uRP:levD_local";
				my $SAMPLE= "0/1:";
				$SAMPLE .= "$plus_Reads[6]:$lSC:$nSC:$num_unmapped_pairs:$local_levD";

				#### remove any white spaces.
				$INFO_1=~s/\s//g;
				$INFO_2=~s/\s//g;

				$BND1_name =~ s/^BND/LEVD/;
				# If left, then the start position is plus_Reads[1]-isize
				my $start_pos=0;
				#Make sure Ref1 and Ref2 bases are different
				if($REF2_base eq $REF1_base){$REF1_base="NA"}
				if($side=~/left/){$start_pos=$plus_Reads[1]-$isize}else{$start_pos=$plus_Reads[1]};		
				print OUT join ("\t", $plus_Reads[0], $start_pos, $BND1_name, $REF2_base, $REF1_base, $QUAL, "PASS", $INFO_1,$FORMAT,$SAMPLE, "\n");
				if ($verbose){print "No Softclipped reads here!\n"}
				return;
			}
		}

		#### Otherwise, look for DRP mate info
	#if($nRP=~/NA/){print "MATE_INFO=$mate_info\tSide=$side\tline=$line\n";}
		my @mate_info_arr = split(/\t/, $mate_info);
		$nRP = $mate_info_arr[3];
		my $mate_chr=$mate_info_arr[0];

			if((! defined $nRP) || ($nRP =~ /na/i) || ($mate_chr =~ /NA/) ){
			#PRINT UNKNOWN
	if ($nRP =~ /na/i){print "Can't find SC reads\n" if ($verbose);return};
	if ($verbose){print "There is an unknown\nNRP=$nRP Mate_CHR=$mate_chr minRP=$minRP\n"}
				$summary .= ":unknown," . $plus_Reads[6];
				$type = "unknown";
				$REF1_base = getSeq($plus_Reads[0], $plus_Reads[1], $INPUT_FASTA);
				$REF2_base = $plus_Reads[3];
				$BND1_name =~ s/^BND/UNKNOWN/;
				$QUAL = 1/($local_levD + 0.001);
				$QUAL = sprintf("%.2f",$QUAL);
				$INFO_1=join(";", "SVTYPE=unknown", "EVENT=unknown", "ISIZE=unknown");
				#Add Sample infomration
				my $FORMAT="GT:sDEL";
				$FORMAT .= ":lSC:nSC:uRP:levD_local";
				my $SAMPLE = "0/1:";
				$SAMPLE .= "$plus_Reads[6]:$lSC:$nSC:$num_unmapped_pairs:$local_levD";
				$SAMPLE=~s/NA/0/g;
				#### remove any white spaces.
				$INFO_1=~s/\s//g;
			       #print join ("\t", $plus_Reads[0], $plus_Reads[1],  $REF2_base, $REF1_base, $QUAL, "PASS", $INFO_1,$FORMAT,$SAMPLE, "\n");

				print OUT join ("\t", $plus_Reads[0], $plus_Reads[1], $BND1_name, $REF1_base, $REF2_base, $QUAL, "PASS", $INFO_1,$FORMAT,$SAMPLE, "\n");
				return;

		}
		#### end if there is no mate info or nRP+uRP<minRP
		if (($nRP<$minRP)&&($unmated_pairs > ($num_unmapped_pairs+$nRP))){
			print "Something failed here\nif (($nRP<$minRP)&&($unmated_pairs > ($num_unmapped_pairs+$nRP)))\n";
		return}

		##################################################################################
		# Find out if mates have nearby soft-clips (to refine the breakpoints)
		##################################################################################
		#Look for evidence of soft-clipping near mate
		my @mate_soft_arr = ();
		my $mate_start = 0;
		my $mate_soft = "";

		@mate_info_arr = split(/\t/, $mate_info);

		#### mate start and end locations.
		my $filename = $right_sc;

		$start = $mate_info_arr[1] - $dist_To_Soft;
		$end = $mate_info_arr[1];

		if ($side =~ /right/) {
			$start = $mate_info_arr[2];
			$end = $mate_info_arr[2] + $dist_To_Soft;

			$filename = $left_sc;
		}

		#### add levenstein distance to Summary
	#print "Calc distal Levd\n";
		$levD->search(rev_comp($plus_Reads[3]), $mate_info_arr[0], $start, $end, $INPUT_FASTA);
		$distl_levD = sprintf("%.2f", $levD->{relative_edit_dist});
	$distl_levD = "NA" if($plus_Reads[3] =~ /NA/);
	#If there is no softclips to string match, then give 0 as quality value
       if ($plus_Reads[3] !~ /NA/){
			$QUAL=1/($distl_levD + 0.001);
		}
		else	{
			$QUAL=0;
		};
	$QUAL=sprintf("%.2f",$QUAL);
	#### looking for softclips to refine break point
	#### if left look in right and vice-versa.
	$cmd = qq{echo -e "$mate_info_arr[0]\t$start\t$end"};
	$cmd .= qq{ | awk -F'\t' 'NF==3' | intersectBed -a stdin -b $filename | head -1};
print "$cmd\n" if $verbose;
	$mate_soft = `$cmd`;

	$mate_soft =~ s/\n//g;
	@mate_soft_arr = split(/\s/, $mate_soft);
my $NO_MATE_SC="";
	if(@mate_soft_arr){
		$mate_chr = $mate_soft_arr[0];
		$mate_start = $mate_soft_arr[1];
                $NO_MATE_SC="APPROXIMATE";

	} else{
		@mate_info_arr = split(/\s/,$mate_info);
		$mate_chr = $mate_info_arr[0];
		$mate_start = $mate_info_arr[1];
	}

	#end if there is no mate info
	return if ($mate_chr eq "");
	#end if there is no mate info and !disable_RP_only
	return if (($lSC =~/NA/)&&(!$disable_RP_only));
	
	
	###########################################################################
	##### Printing output

	#########################################
	# Get DNA info
	#########################################
	#print "PLUS_READS=$plus_Reads[0],$plus_Reads[1]\nMATE=$mate_chr,$mate_start,$INPUT_FASTA\n";
	$REF1_base = getSeq($plus_Reads[0], $plus_Reads[1], $INPUT_FASTA);

	### this is alt ref. for softclip its $plus_Reads[3]
	$REF2_base = getSeq($mate_chr, $mate_start, $INPUT_FASTA);

	#########################################
	# print in VCF format
	#########################################

	#### abs value to account for left and right reads.
	$isize = abs($plus_Reads[1]-$mate_start);
	
	my $event_type=$type;
	$event_type=~ s/,|[0-9]//g;
	$INFO_1=join(";", "$sv_type", "EVENT=$event_type","END=$mate_start", "ISIZE=$isize","MATEID=$BND2_name");
	$INFO_2=join(";", "$sv_type", "EVENT=$event_type","END=$plus_Reads[1]", "ISIZE=$isize","MATEID=$BND1_name");

	#### remove any white spaces.
	#### ask: did you mean to remove space from ends? eg. trim()
	$INFO_1=~s/\s//g;
	$INFO_2=~s/\s//g;

	$FORMAT=$summary;
 	$FORMAT=~ s/,|[0-9]//g;
        $FORMAT .= ":lSC:nSC:uRP:distl_levD";
	if($NO_MATE_SC){$INFO_2 .= ":NO_MATE_SC"}
	my $SAMPLE="0/1:";	
	$SAMPLE .=$summary;
#        if($NO_MATE_SC){$SAMPLE.= ":$NO_MATE_SC"}

	$SAMPLE=~s/[A-Z|,|_]//g;
        my $MATE_SAMPLE=$SAMPLE;
        $SAMPLE .= ":$lSC:$nSC:$num_unmapped_pairs:$distl_levD";
	$MATE_SAMPLE .=":NA:NA:NA:NA";
	$SAMPLE=~s/::/:/g;
	$MATE_SAMPLE=~s/::/:/g;
	$MATE_SAMPLE=~s/NA/0/g;
	$SAMPLE=~s/NA/0/g;
 
	if($type !~ /INV/){
		$ALT_1 = join("","]",$mate_chr,":",$mate_start,"]",$REF1_base);
		$ALT_2 = join("",$REF2_base,"[",$plus_Reads[0],":",$plus_Reads[1],"[");
		#		2      321682 bnd_V  T   ]13:123456]T  6    PASS SVTYPE=BND
		#		13     123456 bnd_U  C   C[2:321682[   6    PASS SVTYPE=BND
	} else {
		$ALT_1 = join("", "]", $plus_Reads[0], ":", $plus_Reads[1], "]", $REF2_base);
		$ALT_2 = join("", $REF1_base, "[", $mate_chr, ":", $mate_start, "[");
	}

	if(($mate_chr) && ($plus_Reads[0])){
		print OUT join ("\t", $plus_Reads[0], $plus_Reads[1], $BND1_name, $REF1_base, $ALT_1, $QUAL,"PASS", $INFO_1, $FORMAT,$SAMPLE,"\n");
		print OUT join ("\t", $mate_chr, $mate_start, $BND2_name, $REF2_base, $ALT_2, $QUAL, "PASS", $INFO_2, $FORMAT,$MATE_SAMPLE,"\n");
	}
}

###############################################################################
###############################################################################
sub get_counts_n_info {
        my ($event, $side, $mapQ, $file, $dist, $target, $upL, $lwL) = @_;

        my $mate_info = "";
        my $cmd = "";

        if ($event =~ /^CTX$/i) {
                #print "CTX side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 32 -F 1536 $file $target};
                        $cmd .= qq{ | perl -ane 'if(\$F[6] ne "="){\$end=\$F[7]+1; print join ("\\t",\$F[6],\$F[7],\$end,"\\n")}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info=`$cmd`;
                } else {
                        $cmd = qq{ samtools view $new_blacklist -q $mapQ -f 16 -F 1536 $file $target};
                        $cmd .= qq{ | perl -ane 'if(\$F[6] ne "="){\$end=\$F[7]+1; print join ("\\t",\$F[6],\$F[7],\$end,"\\n")}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info=`$cmd`;
                }
        } elsif ($event =~ /^DEL$/i) {
                #print "DEL side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ  -f 32 -F 1552 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$9>$upL)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info=`$cmd`;
                } else {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -F 1568 -f 16 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"} {if((\$7 ~ /=/)&&(\$9<-$upL)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist_To_Soft -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}

                        $mate_info=`$cmd`;
                }
        } elsif ($event =~ /^INS$/i) {
                #print "INS side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 32 -F 1552 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$9<$lwL && \$9 > 0 )){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                } else {
                        $cmd = qq {samtools view $new_blacklist -q $mapQ -f 16 -F 1568 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$9>-$lwL && \$9 < 0 )){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                }
        } elsif ($event =~ /^INV$/i) {
                #print "INV side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ  -F 1596 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist_To_Soft -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                } else {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 48 -F 1548 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                }
        } elsif ($event =~ /^TDUP$/i) {
                #print "TDUP side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 32 -F 1552 $file $target};
#			$cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$9>$upL)){end=\$8+1;print \$3,\$8,end}}'};
			$cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$4>\$8)&&(\$9<0)&& (\$9>$upL)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                } else {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 16 -F 1568 $file $target};
#                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$9<-$upL )){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | awk '{OFS="\\t"}{if((\$7 ~ /=/)&&(\$4<\$8)&&(\$9>0)&&(\$9>$upL)){end=\$8+1;print \$3,\$8,end}}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                }
        } elsif ($event =~ /^NOV_INS$/i) {
                #print "NOV_INS side $side\n";
                if ($side =~ /right/i) {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ -f 8 -F 1552 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{end=\$8+1;print \$3,\$8,end}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                } else {
                        $cmd = qq{samtools view $new_blacklist -q $mapQ  -f 24 -F 1536 $file $target};
                        $cmd .= qq{ | awk '{OFS="\\t"}{end=\$8+1;print \$3,\$8,end}'};
                        $cmd .= qq{ | sortBed | mergeBed -d $dist -n | sort -k4nr | head -1};
#if($verbose){print "$cmd\n"}
                        $mate_info = `$cmd`;
                }
        }

        $mate_info=~s/\n//g;
        my @tmp=split(/\t/, $mate_info);

        my $counts = 0;

        if (defined $tmp[3]) {
                $tmp[3] =~ s/\n//g;

                $counts = $tmp[3] if (length($tmp[3]));
        }
        return ({count=>$counts, info=>$mate_info});                                                                                                                                
}
