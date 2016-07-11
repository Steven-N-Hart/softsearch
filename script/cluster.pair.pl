#!/usr/bin/perl                                                                                                                                            
use strict;
use POSIX;

my $usage = "cluster.pair.pl maxdist\n";
my $maxdist = shift or die $usage;

my %count;

while (<STDIN>){
    chomp;
    my ($sample, $chrstart, $start, $chrend, $end) = split /\t/;
    my $nstart = floor ($start/$maxdist);
    my $nend   = floor ($end/$maxdist);
    my $coord = {start=>$start, end=>$end};

    push @{$count{$chrstart}->{$nstart}->{$chrend}->{$nend}->{$sample}}, $coord;
}

print_groups (\%count);

sub print_groups {
    my ($rcount) = @_;
    my %count = %{$rcount};

    foreach my $chrstart (sort {$a<=>$b} keys %count) {
	foreach my $posstart (sort {$a<=>$b} keys %{$count{$chrstart}}) {
	    my %fcoord = %{$count{$chrstart}->{$posstart}};

	    foreach my $chrend (sort {$a<=>$b} keys %fcoord) {
		foreach my $posend (sort {$a<=>$b} keys %{$fcoord{$chrend}}){
		    my @nsamples = sort {$a cmp $b} (keys %{$fcoord{$chrend}->{$posend}});

		    my $cpos = $fcoord{$chrend}->{$posend};

		    my @coords;
		    my $totnum=0;
	    
		    foreach my $sample (@nsamples) {
			my ($num, $avgx, $avgy) = calc_moments(@{$cpos->{$sample}});
			push (@coords, {start=>$avgx, end=>$avgy});
			$totnum+=$num;
		    }

		    my ($num, $avgx, $avgy)  = calc_moments(@coords);
	    
		    print $chrstart."\t".$avgx."\t".$chrend."\t".$avgy ."\t".$num."\t".$totnum."\t" ;
	    
		    print $_."\t" foreach (@nsamples);
		    print "\n";
		}
	    }
	}
    }
}

sub calc_moments {
    my (@pos) = @_;

    my ($num, $sumx, $sumy) = (0,0,0);
    foreach my $cpos (@pos) {
	$num++;
	$sumx+=$cpos->{start};
	$sumy+=$cpos->{end};
    }
    my $avgx = sprintf ("%d", $sumx/$num);
    my $avgy = sprintf ("%d", $sumy/$num);

    return ($num, $avgx, $avgy);
}
