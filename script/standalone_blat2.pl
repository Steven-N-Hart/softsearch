#!/usr/bin/perl -sw
use Getopt::Long;
sub usage(){
    print "
    Usage: <VCF> -g <genome.2bit> -seq|s <seq.fa> -f genome.fa 
	-o out.vcf
	-n contig.names
        -dist   how wide of a window to look for bp [50]\n
	-v	verbose option
        Requires samtools,bedTools, and blat in your path\n;
        ";
    die;
}
#Initialize values
my ($blat,$genome,$tei_bed,$vntr_bed,$out_vcf,$contig_names,$contig,$fasta,$uninformative_contigs,$dist,$verbose,$bedTools,$samtools);
GetOptions ("genome|g=s" => \$genome,
            "o|out:s" => \$out_vcf,
            "names|n:s" => \$contig_names,
            "seq|s=s" => \$contig,
            "f|fasta:s" => \$fasta,
	    "b|bad:s" => \$uninformative_contigs,
            "dist:s" => \$dist,
	    "v" => \$verbose
	    );
#$genome="/data2/bsi/reference/db/hg19.2bit""
#$blat="/projects/bsi/bictools/apps/alignment/blat/34/blat" ;
#TEI.bed=egrep "LINE|SINE|LTR" /data5/bsi/epibreast/m087494.couch/Reference_Data/Annotations/hg19.repeatMasker.bed >TEI.bed
#VNTR_BED=egrep "Satellite|Simple_repeat" /data5/bsi/epibreast/m087494.couch/Reference_Data/Annotations/hg19.repeatMasker.bed > VNTR.bed


$blat=`which blat`;
if (!$blat) {die "Your do not have BLAT in your path\n\n"}
$samtools=`which samtools`;
if (!$samtools) {die "Your do not have samtools in your path\n\n"}
$bedTools=`which sortBed`;
if (!$bedTools) {die "Your do not have bedTools in your path\n\n"}


if (!$dist) {$dist=50}
if (!$out_vcf) {$out_vcf="out.vcf"}
if (!$contig_names) {$contig_names="contig.names"}
if (!$uninformative_contigs) {$uninformative_contigs="uninformative.contigs"}

if ((!$genome)||(!$contig)||(!$fasta)){&usage;die;}


open(VCF,"$ARGV[0]") or die "must specify VCF file\n\n";
open(OUT_VCF,">",$out_vcf) or die "can't open the output VCF\n";
open(CONTIG_LIST,">",$contig_names) or die "can't open the contig names\n";
open(BAD_CONTIG_LIST,">",$uninformative_contigs) or die "can't open the contig names\n";
#print "writing to CONTIG_LIST=$contig_names\n";
while (<VCF>) {
    if($_=~/^#/){
        if ($.==1) {
            print OUT_VCF $_;
            print OUT_VCF "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand to which assembled contig aligned\">\n";
            print OUT_VCF "##INFO=<ID=CONTIG,Number=1,Type=String,Description=\"Name of assembeled contig matching event\">\n";
            print OUT_VCF "##INFO=<ID=MECHANISM,Number=1,Type=String,Description=\"Proposed mechanism of how the event arose\">\n";
            print OUT_VCF "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of insertion\">\n";
            print OUT_VCF "##INFO=<ID=HOM_LEN,Number=1,Type=Integer,Description=\"Length of microhomology\">\n"; 
            next;
        }
    else {
        print OUT_VCF $_;
        next;
        }
    };
    chomp;

    ##look for exact location of BP
    @line=split("\t",$_);
    my($left_chr,$start,$end);

    #Get left position
    $left_chr=$line[0];
    $start=$line[1]-$dist;
    $end=$line[1]+$dist;

    #Get right position
    my ($mate_pos,@mate,$mate_chr,$mate_start,$mate_end);
    $mate_pos=$line[4];
    $mate_pos=~s/[\[|\]|A-Z]//g;
    #print "mate_pos=$mate_pos\n";
    @mate=split(/:/,$mate_pos);
    $mate_chr=$mate[0]; $mate_pos=$mate[1];
    $mate_start=$mate_pos-$dist;$mate_end=$mate_pos+$dist;
    #print "$left_chr:$start-$end\n$mate_chr:$mate_start-$mate_end\n";
    
    #Run through blat
    my ($result1,$result2);
    my $target1=join("",$left_chr,":",$start,"-",$end);
    my $target2=join("",$mate_chr,":",$mate_start,"-",$mate_end);
    #print "target1=$target1\ttarget2=$target2\n";die;
    $result1=get_result($target1);
    $result2=get_result($target2);
   

    my $NOV_INS="";
    #If there is a NOV_INS, then there shouldn't be any output, so trick the results
    if ($_=~/EVENT=NOV_INS/) {
        $mate_start=$start;
        $NOV_INS="true";
        if (!$result1) {$result1=join("\t","0","0","0","0","0","0","0","0","+","UNKNOWN_NODE","0","0",$dist);}
        if (!$result2) {$result2=join("\t","0","0","0","0","0","0","0","0","+","UNKNOWN_NODE","0","0",$dist);}
   }
    
    #Skip over events that aren't supported
    if ((!$result1)||(!$result2)){
	my @tmp1=split("\t",$result1);
	my @tmp2=split("\t",$result2);
	if ($tmp1[9]) {print BAD_CONTIG_LIST "$tmp1[9]\n"}
	if ($tmp2[9]) {print BAD_CONTIG_LIST "$tmp2[9]\n" }
	next;
    }
    #Parse blat results   
    my @result1=split("\t",$result1);
    my @result2=split("\t",$result2);
if($result2[9] ne $result1[9]){print "$result2[9] != $result1[9]\n";next}
    #print "@result1\n@result2\n";die;
    my $pos1=$start+($result1[12]-$result1[11]);
    my $pos2=$mate_start+($result2[12]-$result2[11]);
    #print "$_\n$pos1\t$pos2\n";
    
    ##############################################################
    ### Build Classifier
    
    my ($QSTART1,$QEND1,$QSTART2,$QEND2,$len,$MECHANISM, $INSERTION, $DELETION, $bed_res1,$bed_res2);
    $MECHANISM="UNKNOWN";
    $len="UNKNOWN";
    #Make sure the later event is second
    if ($result1[11] <  $result2[11]){
	$QSTART1=$result1[11];
	$QEND1=$result1[12];
	$QSTART2=$result2[11];
	$QEND2=$result2[12];
    }
    else{
	$QSTART1=$result2[11];
	$QEND1=$result2[12];
	$QSTART2=$result1[11];
	$QEND2=$result1[12];
    }
    #Now calculate the difference between $QEND1 and QSTART2
    if($verbose){print "QEND1=$QEND1\tQSTART2=$QSTART2\n";}
    $len=$QEND1-$QSTART2;
    #Check for TEI
    if($_=~/MECHANISM=TEI/){$MECHANISM="TEI"}
    elsif($_=~/MECHANISM=VNTR/){$MECHANISM="VNTR"}
    else{
        if ($len==0) {$MECHANISM="NHEJ"}
	else{
	    if ($len>0){$INSERTION="true"}
		if ($len<0){$DELETION="true"}
		    if ($INSERTION){
		        if ($len>10) {$MECHANISM="FOSTES"}
		        else{$MECHANISM="NHEJ"}
		    }
		elsif ($DELETION){
		    if ($len>100) {$MECHANISM="NAHR"}
		        elsif ($len > 2){$MECHANISM="altEJ"}
		        else{$MECHANISM="NHEJ"}
	        }
	    }	
	}

    
#if ($verbose){print "@result1";print "@result2";}

    #print out VCF
    #############################################################
    # create temporary variable name
    #############################################################
    srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
    my $random_name=join "", map { ("a".."z")[rand 26] } 1..8;
    my $random_name2=join "", map { ("a".."z")[rand 26] } 1..8;
   
   #Get Ref Base
   my ($ref_base,$alt_base,$tmp_mate_pos);
   $ref_base=getBases($left_chr,$pos1,$fasta);
   $alt_base=getBases($mate_chr,$pos2,$fasta);#print "ALT=$alt_base\n";
   #Substitute the new mate position and base
   $tmp_mate_pos=$line[4];
   $tmp_mate_pos=~s/$mate_pos/$pos2/;
   $tmp_mate_pos=~s/[A-Z]/$alt_base/;
   #split apart the INFO field to adjust the ISIZE and MATEID
   my $NEW_INFO="";
   my @INFO=split(/;/,$line[7]);
   for (my $i=0;$i<@INFO;$i++){
        if ($INFO[$i] =~ /^ISIZE=/){
            my @tmp=split(/=/,$INFO[$i]);
            $NEW_INFO.="ISIZE=";
            my $new_ISZIE=$pos2-$pos1;
            $NEW_INFO.=$new_ISZIE
            }
        elsif($INFO[$i] =~ /^MATE_ID=/){
            $NEW_INFO.=";MATE_ID=".$random_name2 . ";";
        }
        else{
            $NEW_INFO.=$INFO[$i].";";
        }
   }
   #ADD in strand and name
   $NEW_INFO.="STRAND=".$result1[8];
   $NEW_INFO.=";CONTIG=".$result1[9];
   if($MECHANISM!~/TEI|VNTR/){$NEW_INFO.=";MECHANISM=".$MECHANISM;}
    $NEW_INFO.=";HOM_LEN=".$len;
   #don't pring contig nage if its a novel insertion
   if(!$NOV_INS){print CONTIG_LIST "$result1[9]\n";}#else{print "I'm not printing $result1[9]\n";}
    print OUT_VCF "$left_chr\t$pos1\t$random_name\t$ref_base\t$tmp_mate_pos\t1000\tPASS\t$NEW_INFO\t$line[8]\t$line[9]\n";
    #Now go through and fill info in for mate
    #Substitute the new mate position and base
   $tmp_mate_pos=$line[4];
   $tmp_mate_pos=~s/$mate_pos/$pos1/;
   $tmp_mate_pos=~s/[A-Z]/$ref_base/;
   $tmp_mate_pos=~s/$mate_chr/$left_chr/;
    $NEW_INFO="";
    @INFO=split(/;/,$line[7]);
   for (my $i=0;$i<@INFO;$i++){
    if ($INFO[$i] =~ /^ISIZE=/){
            my @tmp=split(/=/,$INFO[$i]);
            $NEW_INFO.="ISIZE=";
            my $new_ISZIE=$pos2-$pos1;
            $NEW_INFO.=$new_ISZIE
            }
        elsif($INFO[$i] =~ /^MATE_ID=/){
            $NEW_INFO.=";MATE_ID=".$random_name.";";
        }
        else{
            $NEW_INFO.=$INFO[$i].";";
        }
   }
    #ADD in strand and name
   $NEW_INFO.="STRAND=".$result2[8];
   $NEW_INFO.=";CONTIG=".$result2[9];
   if ($MECHANISM!~/TEI|VNTR/){$NEW_INFO.=";MECHANISM=".$MECHANISM;}
    $NEW_INFO.=";HOM_LEN=".$len;

   #don't pring contig nage if its a novel insertion
   if(!$NOV_INS){print CONTIG_LIST "$result2[9]\n";} #else{print "I'm not printing $result1[9]\n";}
    print OUT_VCF "$mate_chr\t$pos2\t$random_name2\t$alt_base\t$tmp_mate_pos\t1000\tPASS\t$NEW_INFO\t$line[8]\t$line[9]\n";
	if ($verbose){print  "$mate_chr\t$pos2\t$random_name2\t$alt_base\t$tmp_mate_pos\t1000\tPASS\t$NEW_INFO\t$line[8]\t$line[9]\n";}
}
close VCF;
close OUT_VCF;
close CONTIG_LIST;
close BAD_CONTIG_LIST;
sub get_result{
        my $target=($_[0]);
if($verbose){print "target=$target\n"}#;die;
        my $cmd="blat $genome:$target $contig /dev/stdout -t=dna -q=dna -noHead|egrep -v \"Searched|Loaded\" |head -1";

if ($verbose){print "$cmd\n"}        #print "$cmd\n";die;
        my $result=`$cmd`;
        next if (!$cmd);
        return ($result);
}
sub getBases{
        my ($chr,$pos1,$fasta)=@_;
        my @result=();
        if ($pos1 <0){print "$pos1 is not a valid position (likely caused by circular MT chromosome)\n";$result[1]="NA";};
        @result = `samtools faidx $fasta $chr:$pos1-$pos1`;
        if(!$result[1]){$result[1]="NA"};
        chomp($result[1]);
        return uc($result[1]);
}


