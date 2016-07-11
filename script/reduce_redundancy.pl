open(BUFF,"$ARGV[0]") or die "no input file found\n";
$range="$ARGV[1]";
my %hash;
my %store;
$prev_chr="";
$next=0;
while(<BUFF>)
{
	chomp($_);
	#print "$.\n";
	if($_ !~ m/^#/)
	{
		@array=split("\t",$_);
		$chr=$array[0];
		$pos=$array[1];
		$value=$array[@array-1];
		if($prev_chr ne $chr )
		{
			if($prev_chr ne "")
			{
				foreach $key (sort {$hash{$b} <=> $hash{$a} } keys %hash)
                        	{
                                	print "$store{$key}\n";
                                	last;
                        	}

			}
			$next = $pos+$range;
			undef(%hash);
			undef(%store);
		}
		if($next< $pos)
		{	
			foreach $key (sort {$hash{$b} <=> $hash{$a} } keys %hash)
			{
     				print "$store{$key}\n";
				last;
			}
			$next = $pos+$range;
			undef(%hash);
			undef(%store);
			
		}	
		if($value eq "NA")
                {
                      $hash{$chr." ".$pos." ".$.}=0;
                }
                else
                {
                       $hash{$chr." ".$pos." ".$.}=$value;
               	}
                $store{$chr." ".$pos." ".$.}=$_;
	}
	else
	{
		print $_."\n";
	}
	$prev_chr = $chr;
}
foreach $key (sort {$hash{$b} <=> $hash{$a} } keys %hash)
{
       print "$store{$key}\n";
       last;
}

