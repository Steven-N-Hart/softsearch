package LevD;

use lib "/data2/bsi/reference/softsearch/lib/perl5";
use strict;
use warnings;
use Data::Dumper;
use String::Approx 'adist';
use String::Approx 'adistr';
use String::Approx 'aindex';

my $WINDOW_SIZE = 100;

sub new {
	my ($class, $file) = @_;
    my $self = {};

 	bless($self,$class);
	$self->init();

	return $self;
}

sub init {
	my ($self) = @_;

	#### default values.
	$self->{index} = 0;
	$self->{relative_edit_dist} = 0;
	$self->{edit_dist} = 0;
}

sub search {
	my ($self, $clip, $chr, $start, $stop, $ref) = @_;

	if (! -s $ref) {
		die "ERROR: Reference file $ref now found\n";
	}

	#### extact seq from reference file.
	my $target = $chr .":". $start ."-". $stop;
	my $cmd = "samtools faidx $ref $target";

	my @output = $self->_run_system_cmd($cmd);

	#### depending on ref file format seq could be on multiple lines
	#### concatinate all except for the header in one line.
	#### e.g:
	#### >chr1:8222999-8223099
	#### GGTGCAATCATAGCTCACTAAGCTTCAACCTCAAGAGATCCTCCCACCTCAGCCTCCCAG
	#### GTAGCTGGGACTACAGGCAAATGCCATGACACCTAGCTAAT
	my $seq = join("", @output[1..$#output]);

	#### remove new line character
	$seq =~ s/\n//g;

	#### find number of mismatches and start index
	#### of clip to be searched against target seq.
	$self->{relative_edit_dist} = adistr($clip, $seq);
	$self->{edit_dist} = adist($clip, $seq);
	$self->{index} = aindex($clip, $seq);
}

sub _run_system_cmd {
	my ($self, $cmd) = @_;
	my @cmd_output;

	eval {
		@cmd_output = qx{$cmd 2>&1};
		if ( ($? << 8) != 0 ) {
			die "@cmd_output";
		}
	};
	if ($@) {
		die "Error executing command $cmd: $@";
	}

	return @cmd_output;
}

1;
