#!/usr/bin/perl

=head1 NAME
   install.pl

=head1 SYNOPSIS
    USAGE: install.pl --prefix=/location/of/install/dir

=head1 OPTIONS

B<--prefix, -p>
	Required. Prefix location where package will be installed.

B<--perl_exec, -e>
	Optional.  If perl exec is other than /usr/bin/perl please specify location of perl install

B<--help,-h>


=head1  DESCRIPTION
	Install package

=head1  INPUT

=head1  OUTPUT


=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./install.pl --prefix=/prefix

=cut

use strict;
use warnings;
use Cwd;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'prefix|p=s',
						  'perl_exec|e=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#############################################################################
#### make sure everything passed was peachy
&check_parameters(\%options);

#### print time now.
timestamp();

my $this = {};
my $progress = {};
my $cmd = "";

#### get current working dir
$this->{source} = getcwd();

$progress = getProgress();

#### make logs dir
$cmd = "mkdir -p $options{prefix}/logs";
execute_cmd($cmd);

#### installling libraries required for successfull run
install_libraries();

#### unpack binary dir containing all binary to be installed
#### which are required for successfull run
print STDERR "\n\nInstalling binaries...\n";

#### install each package in binary folder.
my @packages = qw(stringApprox levD);

foreach my $tool (@packages) {
	if ((exists $progress->{$tool}) && ($progress->{$tool})){
		print STDERR "\t$tool already installed. Skipping...\n";
	} else {
		print STDERR "\tInstalling $tool...\n";

		#### unpack and install each tool
		eval("install_${tool}()");
	}
}

#### copy source code and update paths for perl and libs
install_source();

#### completion message
print "\n\n\tSoftSearch installation complete.  Use following command to initiate a test run\n";
print "\n\tperl $options{prefix}/src/SoftSearch.pl -f {GENOME} -b {BAM_FILE}\n\n";

#### print time now
timestamp();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(prefix);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	$options{'perl_exec'} = "/usr/bin/perl" unless($options{'perl_exec'});
}

#############################################################################
sub getProgress {
	my $hash = {};
	my @sofar;

	#### if file exists get progress so far.
	if (-s "$options{prefix}/progress.txt") {
		open(FHD, "<", "$options{prefix}/progress.txt") or die "Could not open file to read $options{prefix}/progress.txt";
		while(<FHD>){
			chomp $_;
			push @sofar, $_;
		}
		close(FHD);

		map { $hash->{$1} = $2 if( /([^=]+)\s*=\s*([^=]+)/ ) } @sofar;
	}

	#### return hash
	return $hash;
}

#############################################################################
sub setProgress {
	my $hash = shift;

	open(OUT, ">", "$options{prefix}/progress.txt") or die "Could not open file to write $options{prefix}/progress.txt";

	foreach my $key (keys %{$hash}){
		print OUT $key."=".$hash->{$key}."\n";
	}

	close(OUT);
}

#############################################################################
sub install_libraries {
	if ((exists $progress->{libraries}) && ($progress->{libraries})){
		print STDERR "\tLibraries already installed. Skipping...\n";
		return;
	}

	print STDERR "\n\nInstalling libraries...\n\n";
	chdir($this->{source});

	$cmd = "cp -r $this->{source}/library $options{prefix}/lib";
	execute_cmd($cmd);

	$progress->{libraries} = 1;
	setProgress($progress);
}

#############################################################################
sub install_stringApprox {
	#### check and install dir
	my $dir = "$options{prefix}/lib";
	my $cmd = "";

	$cmd = "mkdir -p $dir";
	execute_cmd($cmd);

	$cmd = "tar -zxvf $this->{source}/binary/String-Approx-3.27.tar.gz -C $this->{source}/binary";
	execute_cmd($cmd);

	chdir("$this->{source}/binary/String-Approx-3.27");
	$cmd = "perl Makefile.PL INSTALL_BASE=$options{prefix}";
	$cmd .= " 1>$options{prefix}/logs/StringApprox.out";
	$cmd .= " 2>$options{prefix}/logs/StringApprox.err";
	execute_cmd($cmd);

	$cmd = "make && make install";
	$cmd .= " 1>>$options{prefix}/logs/StringApprox.out";
	$cmd .= " 2>>$options{prefix}/logs/StringApprox.err";
	execute_cmd($cmd);

	$cmd = "make install";
	$cmd .= " 1>>$options{prefix}/logs/StringApprox.out";
	$cmd .= " 2>>$options{prefix}/logs/StringApprox.err";
	execute_cmd($cmd);


	chdir("$this->{source}/binary");
	$cmd = "rm -rf $this->{source}/binary/String-Approx-3.27";
	execute_cmd($cmd);

	$progress->{stringApprox} = 1;
	setProgress($progress);
}

#############################################################################
sub install_levD {
	#### check and install dir
	my $dir = "$options{prefix}/lib";
	my $cmd = "";

	$cmd = "mkdir -p $dir";
	execute_cmd($cmd);

	$cmd = "tar -zxvf $this->{source}/binary/Text-LevenshteinXS-0.03.tar.gz -C $this->{source}/binary";
	execute_cmd($cmd);

	chdir("$this->{source}/binary/Text-LevenshteinXS-0.03");
	$cmd = "perl Makefile.PL INSTALL_BASE=$options{prefix}";
	$cmd .= " 1>$options{prefix}/logs/levD.out";
	$cmd .= " 2>$options{prefix}/logs/levD.err";
	execute_cmd($cmd);

	$cmd = "make";
	$cmd .= " 1>>$options{prefix}/logs/levD.out";
	$cmd .= " 2>>$options{prefix}/logs/levD.err";
	execute_cmd($cmd);

	$cmd .= "make install";
	$cmd .= " 1>>$options{prefix}/logs/levD.out";
	$cmd .= " 2>>$options{prefix}/logs/levD.err";
	execute_cmd($cmd);

	chdir("$this->{source}/binary");
	$cmd = "rm -rf $this->{source}/binary/Text-LevenshteinXS-0.03";
	execute_cmd($cmd);

	$progress->{levD} = 1;
	setProgress($progress);
}

#############################################################################
sub install_source {
	if ((exists $progress->{source}) && ($progress->{source})){
		print STDERR "\tSource already installed. Skipping...\n";
		return;
	}

	print STDERR "\n\nInstalling source...\n\n";

	#### create dir to store source code
	$cmd = "mkdir -p $options{prefix}/src";
	execute_cmd($cmd);

	$cmd = "cp -r $this->{source}/script/* $options{prefix}/src/.";
	execute_cmd($cmd);

	#### make sure all scripts are executable
	$cmd = "chmod -R +x $options{prefix}/src";
	execute_cmd($cmd);

	#### replace /usr/local/biotools/perl/5.10.0/bin/perl with perl_exec
	$options{perl_exec} =~ s/\//\\\//g;
	$cmd = "find $options{prefix}/src -name \"*.pl\" -print";
	$cmd .= " -exec sed -i 's/#!\\/usr\\/local\\/biotools\\/perl\\/5.10.0\\/bin\\/perl/#!$options{perl_exec}/' {} \\;";
	execute_cmd($cmd);

	#### check if perl exec location is other than /usr/bin/perl
	if ($options{perl_exec} !~ /^\/usr\/bin\/perl$/) {
		$cmd = "find $options{prefix}/src -name \"*.pl\" -print";
		$cmd .= " -exec sed -i 's/#!\\/usr\\/bin\\/perl/#!$options{perl_exec}/' {} \\;";
		execute_cmd($cmd);
	}

	#### replace library references to local install
	my $lib = "$options{prefix}/lib";
	$lib =~ s/\//\\\//g;

        $cmd = "find $options{prefix}/src -name \"*.pl\" -print";
        $cmd .= " -exec sed -i 's/\\/data2\\/bsi\\/reference\\/softsearch\\/lib/$lib/' {} \\;";
        execute_cmd($cmd);
        
        $cmd = "find $options{prefix}/lib -name \"LevD.pm\" -print";
        $cmd .= " -exec sed -i 's/\\/data2\\/bsi\\/reference\\/softsearch\\/lib/$lib/' {} \\;";
        execute_cmd($cmd);

	$progress->{source} = 1;
	setProgress($progress);
}

#############################################################################
sub execute_cmd {
	my $cmd = shift;

	system($cmd);

	#while (( $? >> 8 ) != 0 ){
	#	print STDERR "ERROR: Following command failed to execute. Exiting execution of workflow\n$cmd\n";
	#	exit(-1);
	#}
}

#############################################################################
sub timestamp {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime."\n";
}
