#####################################################################################################################################################
#Purpose: To parse blat psl file
#Date: 07-30-2013
#####################################################################################################################################################
use Getopt::Long;
use Cwd;
#reading input arguments
&Getopt::Long::GetOptions(
'b|BLAT_OUT=s'=> \$blat_out,
'temp:s'=>\$dirtemp,
'f|FASTA=s'=>\$infast,
);
$blat_out =~ s/\s|\t|\r|\n//g;
$dirtemp =~ s/\s|\t|\r|\n//g;
$infast =~ s/\s|\t|\r|\n//g;
$samtools=`which samtools`;
$samtools =~ s/\s|\t|\r|\n//g;

if($blat_out eq "" || $infast eq "" )
{
	die "Try: perl blat_parse.pl -b <PSL FILE> -f <Contigs.fa> 
	-temp	temporary file directory
	\n";
}   
if (!(-e $samtools))
{
	die "samtools must be in your path\n";
}

if (!(-e $infast))
{
	die "input fasta file doesn't exit\n";
}
unless(-d $dirtemp)
{
    #system("mkdir -p $dirtemp");
    $dirtemp= getcwd;
}	
#opening the blat output file
open(BUFF,$blat_out) or die "no file found $blat_out\n";
open(WRBUFF,">$dirtemp/Temp_out.txt") or  die "not able to write the file \n";
#parsing throught he file
while(<BUFF>)
{
	if($_ =~ m/^\d/)
	{
		print WRBUFF $_;	
	}
	else
	{
		print "ignoring headers $.\n";
	}
}	
close(WRBUFF);
system("sort -k10,10 -k18,18n $dirtemp/Temp_out.txt > $dirtemp/Temp_out1.txt");
system("mv  $dirtemp/Temp_out1.txt $dirtemp/Temp_out.txt");
open(BUFF,"$dirtemp/Temp_out.txt") or die "no file found Temp_out.txt\n";
open(WRBUFF,">$dirtemp/File1_out.txt") or  die "not able to write the file \n";
close(WRBUFF);

$prev_contig_name="";
my @temp;
#parsing throught he file
while(<BUFF>)
{
	
		chomp($_);
		split "\t";
		if($_[9] ne $prev_contig_name)
		{
			if($prev_contig_name ne "")
			{
				#print @temp."\n";
				#print @temp."\n";
				&processing(@temp);
			}
			undef(@temp);
			push(@temp,$_);		
		}
		else
		{
			push(@temp,$_);
		}	
		$prev_contig_name=$_[9];	
	
	
}	
#processing last record
&processing(@temp);
#print @temp."\n";
close(BUFF);




##################SUBROUTINES######################
#actual processing of each record in the temp array(same query name objects)

sub processing {
	open(WRBUFF,">>$dirtemp/File1_out.txt") or  die "not able to write the file \n";
        open(BAD_CONTIG,">>$dirtemp/bad_contig.out.txt") or  die "not able to write the file \n";

	@temp = @_;
	#if number of hits for a contig is one
	if(@temp == 1)
	{
			$i=0;
			#define blocksizes array
			@row=split("\t",$temp[$i]);
			$row[18] =~ s/,$//g;
			@blockSizes=split(',',$row[18]);
			#defining var
			$qSize=$row[10];
			$qStart=$row[11];
			$qStop=$row[12];
			$tstart=$row[15];
			$tstop=$row[16];
			$Strand=$row[8];
			$coverage = $row[9];
			$coverage =~ s/\w+_//g;
			#calculate match val
			if(($qSize-($qStop-$qStart)) ==0)
			{ 	
				$flag=1;
				#these ara non informative
				if (@blockSizes ==1)
				{
					print "ignoring one of the event $row[9] $i as the event is non informative \n";
					print BAD_CONTIG "$row[9]\n";
				}
				#Ignoring when number of blocks are more than two
				if(@blockSizes > 2)
				{
					print "ignoring event $row[9] $. AS BLOCK SIZE is greater than 2\n";	
				}
				#if number of blocks is equal to 2
				if(@blockSizes == 2)
				{
					$temp1=$tstart+$blockSizes[0]+1;
					$temp2=$tstop-$blockSizes[1]-1;
						
					print  WRBUFF "$row[9]\t$row[13]\t$temp1\t$Strand\t$row[13]\t$temp2\t$Strand\t$coverage\n";
				}
				$i=@temp;
			}
			#later part missing
			elsif($qStart ==0)
			{	
				$temp1=$tstart+$blockSizes[0]+1;
				$infast_chr=$infast;
				$infast_chr=~ s/\.fa//g;
				$infast_chr_start=$qStop+1;
				$infast_chr_stop=$qSize;
				$sys="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
				
				$sys = `$sys`;
				chomp($sys);
				@sys=split("\n",$sys);
				$INSERTION="";
				for($i=1;$i<@sys;$i++)
				{
					$INSERTION=$INSERTION.$sys[$i];
				}
				$INSERTION_LENGTH=length($INSERTION);
				$temp1=$tstart+$blockSizes[0]+1;
				print  WRBUFF "$row[9]\t$row[13]\t$temp1\t$Strand\tUNKNOWN\tUNKNOWN\t$Strand\t$coverage\t$INSERTION\t$INSERTION_LENGTH\n";
				
			}
			#intial part missing
			elsif($qStop == $qSize)
			{
				$temp1=$tstart;
				$infast_chr=$infast;
				$infast_chr=~ s/\.fa//g;
				$infast_chr_start=0;
				$infast_chr_stop=$qStart;
				$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
				#die "$sys\n";
				$sys = `$sys`;
				#die "$sys\n";
				chomp($sys);
				@sys=split("\n",$sys);
				$INSERTION="";
				for( $i=1;$i<@sys;$i++)
				{
						$INSERTION=$INSERTION.$sys[$i];
				}
				$INSERTION_LENGTH=length($INSERTION);
				$temp1=$tstart+1;
				print  WRBUFF "$row[9]\tUNKNOWN\tUNKNOWN\t$Strand\t$row[13]\t$temp1\t$Strand\t$coverage\n";
				
			}
			else
			{
				print "ignoring one of the event $row[9] $i as the event is non informative \n";
			}
		
	}
	#if number of hits for a contig is greater than one
	else
	{
		#this flag is used to see if perfect hit not found (match val =0)
		$flag1 = 0;
		for(my $i=0;$i<@temp;$i++)
		{
			
			#define blocksizes array
			@row=split("\t",$temp[$i]);
			$row[18] =~ s/,$//g;
			@blockSizes=split(',',$row[18]);
			#defining var
			$qSize=$row[10];
			$qStart=$row[11];
			$qStop=$row[12];
			$tstart=$row[15];
			$tstop=$row[16];
			$Strand=$row[8];
			$coverage = $row[9];
			$coverage =~ s/\w+_//g;
			#calculate match val
			if(($qSize-($qStop-$qStart)) ==0)
			{ 	
				$flag1=1;
				#these ara non informative
				if (@blockSizes ==1)
				{
					print "ignoring one of the event $row[9] $i as the event is non informative \n";
					print BAD_CONTIG "$row[9]\n";
				}
				#Ignoring when number of blocks are more than two
				if(@blockSizes > 2)
				{
					print "ignoring event $row[9] $. AS BLOCK SIZE is greater than 2\n";	
				}
				if(@blockSizes == 2)
				{
					$temp1=$tstart+$blockSizes[0]+1;
					$temp2=$tstop-$blockSizes[1]-1;
						
					print  WRBUFF "$row[9]\t$row[13]\t$temp1\t$Strand\t$row[13]\t$temp2\t$Strand\t$coverage\n";
				}
				$i=@temp;
			}
		}
		#as flag value not changed proceed to see next step
		if($flag1 == 0)
		{
			undef(@initial);
			my @initial;
			for(my $i=0;$i<@temp;$i++)
			{
				@row=split("\t",$temp[$i]);
				#print "@row\n";
				unshift(@initial,[@row]);
			}
			#sortin the hits according to qstart & qend
			@initial = sort {$a->[11] <=> $b->[11] || $b->[12] <=> $a->[12]} @initial;
			#print "$row[9]\t@initial\n";
			#if($row[9]  eq "NODE_5_length_149_cov_12.395973")
			#{
			#	for($i=0;$i<@initial;$i++)
			#	{
			#		print "@{$initial[$i]}\n";
			#	}
			#}
			$start = "";
			$stop = "";
			$start_len=0;
			$stop_len=0;
			#this super flag is used to skip processing of remaining uncessary hits
			$super_flag = 0;
			for($i=0;$i<@initial && $super_flag == 0;$i++)
			{
				$flag = 0;
				#print "@{$initial[$i]}\n";
				$initial[$i][18] =~ s/,$//g;
				@blockSizes1=split(',',$initial[$i][18]);
				#defining var
				$qSize1=$initial[$i][10];
				$qStart1=$initial[$i][11];
				$qStop1=$initial[$i][12];
				$tstart1=$initial[$i][15];
				$tstop1=$initial[$i][16];
				$Strand1=$initial[$i][8];
				$Chr1 = $initial[$i][13];
				$coverage1 = $initial[$i][9];
				$coverage1 =~ s/\w+_//g;
				#die "$qSize1\t$qStart1\t$qStop1\t$tstart1\t$tstop1\t$Strand1\t$Chr1\t$coverage1\n";
				#if a hit qstart = 0 then set flag =1 
				if($qStart1 == 0)
				{
					$flag =1;
				}
				#if a hit qstop = 0 then set flag =2 
				if($qStop1 == $qSize1)
				{
					$flag =2;
				}
				#if($row[9]  eq "NODE_5_length_149_cov_12.395973")
				#{
				#	print "$flag \n";
				#}
				if(@blockSizes1 == 1)
				{
					if($flag == 1 )
					{
						for($j=0;$j<@initial;$j++)
						{
							#both hits should not be the same 
							if($i != $j)
							{
								#print "@{$initial[$i]}\n";
								$initial[$j][18] =~ s/,$//g;
								@blockSizes2=split(',',$initial[$j][18]);
								#defining var
								$qSize2=$initial[$j][10];
								$qStart2=$initial[$j][11];
								$qStop2=$initial[$j][12];
								$tstart2=$initial[$j][15];
								$tstop2=$initial[$j][16];
								$Strand2=$initial[$j][8];
								$coverage2 = $initial[$j][9];
								$Chr2 = $initial[$j][13];
								$coverage2 =~ s/\w+_//g;
								#making sure both hits are not over lapping
								if($qStart2 > $qStart1)
								{	#allowing +-2 bases as the this hit is immediate next continous hit
									if($qStop1 >= $qStart2 -2  &&  $qStop1 <= $qStart2 +2  )
									{
										#perfect match
										if($qStop2 == $qSize2)
										{
											if($Strand1 eq "+")
											{
												$tmp1 = $tstart1+$blockSizes1[0]+1;
												$tmp2 = $tstart2+$blockSizes2[0];
												print WRBUFF "$initial[$i][9]\t$Chr1\t$tmp1\t$Strand1\t$Chr2\t$tmp2\t$Strand2\t$coverage1\n";
											}
											else
											{
												$tmp1 = $tstart1+1;
												$tmp2 = $tstart2+1;
												print WRBUFF "$initial[$i][9]\t$Chr1\t$tmp1\t$Strand1\t$Chr2\t$tmp2\t$Strand2\t$coverage1\n";
											
											}
											$super_flag = 1;
											$j = @initial+1;	
										}
										#some part is missing after the second hit
										else
										{
											$tmp1 = $tstart1+$blockSizes1[0];
											$tmp2 = $tstart2+$blockSizes2[0];
											$INSERTION="";
											$infast_chr=$infast;
											$infast_chr=~ s/\.fa//g;
											$infast_chr_start=$qStop1+1;
											$infast_chr_stop=$qStart2-1;
											$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
											#die "$sys\n";
											$sys = `$sys`;
											#die "$sys\n";
											chomp($sys);
											@sys=split("\n",$sys);
											for( $i=1;$i<@sys;$i++)
											{
												$INSERTION=$INSERTION.$sys[$i];
											}
											$INSERTION_LENGTH=length($INSERTION);
											print WRBUFF "$initial[$i][9]\t$Chr1\t$tmp1\t$Strand1\t$Chr2\t$tmp2\t$Strand2\t$coverage1\n";
											$super_flag = 1;
											$j = @initial+1;	 
										}
										
									}
									#if there are some insertion between two hits
									elsif($qStop2 == $qSize2)
									{
										$tmp1 = $tstart1+$blockSizes1[0];
										$tmp2 = $tstart2+$blockSizes2[0];
										$INSERTION="";
										$infast_chr=$infast;
										$infast_chr=~ s/\.fa//g;
										$infast_chr_start=$qStop2+1;
										$infast_chr_stop=$qSize;
										$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
										#die "$sys\n";
										$sys = `$sys`;
										#die "$sys\n";
										chomp($sys);
										@sys=split("\n",$sys);
										for( $i=1;$i<@sys;$i++)
										{
											$INSERTION=$INSERTION.$sys[$i];
										}
										$INSERTION_LENGTH=length($INSERTION);
										print WRBUFF "$initial[$i][9]\t$Chr1\t$tmp1\t$Strand1\t$Chr2\t$tmp2\t$Strand2\t$coverage1\n";
										$super_flag = 1;
										$j = @initial+1;	
									}
												
								}
									
							}	
						}
						#if none worked with other reads then only process that read
						if($j == @initial)
						{
							#die "success\n";
							$temp1=$tstart1+$blockSizes1[0]+1;
							#print  WRBUFF "$Chr1\t$temp1\t$Strand1\tUNKNOWN\tUNKNOWN\t$Strand\t$coverage\n";
							$infast_chr=$infast;
							$infast_chr=~ s/\.fa//g;
							$infast_chr_start=$qStop1+1;
							$infast_chr_stop=$qSize1;
							$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
							#die "$sys\n";
							$sys = `$sys`;
							#die "$sys\n";
							chomp($sys);
							@sys=split("\n",$sys);
							$INSERTION="";
							for( $i=1;$i<@sys;$i++)
							{
								$INSERTION=$INSERTION.$sys[$i];
							}
							$INSERTION_LENGTH=length($INSERTION);
							print WRBUFF "$initial[$i][9]\t$Chr1\t$temp1\t$Strand1\tUNKNOWN\tUNKNOWN\t$Strand1\t$coverage1\n";
							$super_flag = 1;
						}	
					}
					#if query end is matched to query size
					elsif($flag == 2)
					{
						#going through other hits
						for($j=0;$j<@initial;$j++)
						{
							#hits should not be same
							if($i != $j && $qStop2)
							{
								#print "@{$initial[$i]}\n";
								$initial[$j][18] =~ s/,$//g;
								@blockSizes2=split(',',$initial[$j][18]);
								#defining var
								$qSize2=$initial[$j][10];
								$qStart2=$initial[$j][11];
								$qStop2=$initial[$j][12];
								$tstart2=$initial[$j][15];
								$tstop2=$initial[$j][16];
								$Strand2=$initial[$j][8];
								$coverage2 = $initial[$j][9];
								$Chr2 = $initial[$j][13];
								$coverage2 =~ s/\w+_//g;
								#if 
								if($qStop2 < $qStop1)
								{
									if($qStart1 >= $qStop2 -2  &&  $qStart1 <= $qStop2 +2  )
									{
										#die "$qStart1 <= $qStop2 \n";
										$tmp1 = $tstart1+$blockSizes1[0];
										$tmp2 = $tstart2+$blockSizes2[0];
										$INSERTION="";
										$infast_chr=$infast;
										$infast_chr=~ s/\.fa//g;
										$infast_chr_start=0;
										$infast_chr_stop=$qStart1-1;
										$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
										#die "test $sys\n";
										$sys = `$sys`;
										#die "$sys\n";
										chomp($sys);
										@sys=split("\n",$sys);
										for( $i=1;$i<@sys;$i++)
										{
											$INSERTION=$INSERTION.$sys[$i];
										}
										$INSERTION_LENGTH=length($INSERTION);
										print WRBUFF "$initial[$i][9]\t$Chr2\t$tmp2\t$Strand2\t$Chr1\t$tmp1\t$Strand1\t$coverage1\n";
										$super_flag = 1;
										$j = @initial+1;
										
									}
									
								}	
							}
						}
						if($j == @initial)
						{
							$infast_chr=$infast;
							$infast_chr=~ s/\.fa//g;
							$infast_chr_start=0;
							$infast_chr_stop=$qStart1;
							$sys ="$samtools faidx $infast $infast_chr:$infast_chr_start-$infast_chr_stop";
							#die "test $sys\n";
							$sys = `$sys`;
							#die "$sys\n";
							chomp($sys);
							@sys=split("\n",$sys);
							$INSERTION="";
							for( $i=1;$i<@sys;$i++)
							{
								$INSERTION=$INSERTION.$sys[$i];							
							}
							$INSERTION_LENGTH=length($INSERTION);
							$tmp = $tstart1+1;
							print WRBUFF "$initial[$i][9]\tUNKNOWN\tUNKNOWN\t$Strand1\t$Chr1\t$tmp\t$Strand1\t$coverage1\n";
							$super_flag = 1;
						}	
					}
				}
				elsif(@blockSizes == 2)
				{
					$temp1=$tstart1+$blockSizes[0]+1;
					$temp2=$tstop1-$blockSizes[1]-1;
					print  WRBUFF "$initial[$i][9]\t$Chr1\t$temp1\t$Strand1\t$Chr1\t$temp2\t$Strand1\t$coverage1\n";
				
				}		
			}
		}
		
	}
	close(WRBUFF);
	
	undef(@temp);
}
 
