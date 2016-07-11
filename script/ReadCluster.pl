#!/usr/bin/perl

=head1 NAME
   ReadCluster.pl

=head1 SYNOPSIS

    USAGE: ReadCluster.pl --input input_sam_file --output output_prefix [--threshold 10000 --minClusterSize 4]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--output,-o>
   output prefix

B<--window, -w>
    Window size

B<--minClusterSize, -m>
	Min size of cluster

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   ReadCluster.pl --input=filename.sam --window=10000 --output=PREFIX

=cut

use strict;
use warnings;
use Data::Dumper;
use DBI;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
                          'window|w=s',
						  'minClusterSize|m=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my $r1_start = 0;
my $r2_start = 0;
my $r1_end = $r1_start + $options{window};
my $r2_end = $r2_start + $options{window};
my $r1_chr = "";
my $r2_chr = "";

my @cluster = ();

open (FHD, "<", $options{input}) or die "Cound not open file $options{input}\n";
open (INTRA, ">", $options{output} . ".intra.sam") or die "Cound not open file $options{output}.intra.sam\n";
open (INTER, ">", $options{output} . ".inter.sam") or die "Cound not open file $options{output}.inter.sam\n";

while (<FHD>){
	chomp $_;

	#skip processing lines starting with @ just print to output file.
	if ($_ =~ /^@/){
		print INTRA $_."\n";
		print INTER $_."\n";
		next;
	}
#print "$_\n";
	check_sequence($_);
}

close(FHD);
close(INTRA);
close(INTER);

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{window}) { $options{window} = 10000; }
	unless($options{minClusterSize}) { $options{minClusterSize} = 4; }
}

#############################################################################
sub check_sequence {
	my $line = shift;

	my @data = split(/\t/, $line);

	## check if mates are within the window.
	if ((inWindow($data[3], 1)) && (inWindow($data[7], 2)) &&
		($r1_chr =~ /$data[2]/) && ($r2_chr =~ /$data[6]/)) {

		## if minClusterSize is reached output
		if (scalar(@cluster) >= $options{minClusterSize}) {

			## if chr are the same then print intra-chr else inter-chr
			if ($data[6] =~ /=/) {
				print INTRA $line."\n";
			} else {
				print INTER $line."\n";
			}
		} else {
			push @cluster, $line;
		}
	} else {

		if (scalar(@cluster) >= $options{minClusterSize}) {
			dumpCluster(@cluster);
		}

		@cluster = ();
		$r1_start = $data[3];
		$r2_start = $data[7];
		$r1_end = $r1_start + $options{window};
		$r2_end = $r2_start + $options{window};
		$r1_chr = $data[2];
		$r2_chr = $data[6];
	}
}

#############################################################################
sub inWindow {
	my $coord = shift;
	my $read = shift;

	my $start = 0;
	my $end = 0;

	if ($read == 1) {
		$start = $r1_start;
		$end = $r1_end;
	} else {
		$start = $r2_start;
		$end = $r2_end;
	}

	if (($coord > $start) && ($coord < $end)){
		return 1;
	} else { return 0; }
}

#############################################################################
sub dumpCluster {
	my @cluster = shift;

	foreach (@cluster){
		my @data = split(/\t/, $_);

		if ($data[6] =~ /=/) {
			print INTRA $_."\n";
		} else {
			print INTER $_."\n";
		}
	}
}
