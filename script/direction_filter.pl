use Getopt::Long;
my ($v);

GetOptions ("v|verbose"  => \$v);   # flag



open (FILE,"$ARGV[0]") or die "Cant find file\n\n";
my $dist=0;
my $pos=0;
my @max=0;
my @events=0;

while(<FILE>){
	$dist=0;
	@first=split(/\s+/,$_);
	$numEvents=($_=~tr/\|//)+1;
	$dist=$first[1]-$pos;
	push(@max,$_);
	push(@events,$numEvents);
#print "STARTING_POS=$pos\n";
	if(($dist<500)||eof()){
		until (($dist>500)||eof()){
			$newline=<FILE>;
			@second=split(/\s+/,$newline);
			$numEvents=($newline=~tr/\|//)+1;
			push(@max,$newline);
			push(@events,$numEvents);
			if($v){print "DIST=$dist\nSEC1=@second[1] POS1=$pos;\n";}
			my $tmp=$pos;
			$pos=@second[1];
			$dist=@second[1]-$tmp;
		}
	}
if ($v){print "Corrected dist= $dist\n" if ($v)};
	#Get the last values since they don't count
	$NL=pop(@max);
	$NE=pop(@events);
	my $idxMax = 0;
	#Get the index of the largest value in array
	if ($v){print "Picking from events:\n"};
	$events[$idxMax] > $events[$_] or $idxMax = $_ for 1 .. $#events;

	my $val=@max[$idxMax];
	print "$val" unless ($val=~/^0$/) ;
	
	
	@max=$NL;
	@events=$NE;		
	my @tmp=split(/\s+/,$NL);
	$pos=$tmp[1];
}

close FILE;

