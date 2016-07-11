#####################################################################################################################################################
#Purpose: To perform blat and organize blat
#Date: 07-19-2013
#####################################################################################################################################################
use Getopt::Long;
#reading input arguments
&Getopt::Long::GetOptions(
'BLAT_PATH=s'=> \$blatpath,
'REF_FILE=s'=> \$reffile,
'INPUT_FILE=s' => \$inputfile,
'OUTPUT_FILE=s' => \$outputfile,
'MIN_SCORE=s'=> \$minScore,
'MIN_IDENTITY=s'=> \$minidentity,
'BLAT_PORT=s'=>\$blatport
);
$blatpath =~ s/\s|\t|\r|\n//g;
$reffile=~ s/\s|\t|\r|\n//g;
$inputfile=~ s/\s|\t|\r|\n//g;
$outputfile=~ s/\s|\t|\r|\n//g;
$minScore=~ s/\s|\t|\r|\n//g;
$minidentity=~ s/\s|\t|\r|\n//g;
$blatport=~ s/\s|\t|\r|\n//g;
#input arguments

#checking for missing arguments
if($blatport eq "" || $blatpath eq "" || $reffile eq "" || $inputfile eq "" || $outputfile eq "" || $minScore eq "" || $minidentity eq "")
{
	die "missing arguments\n USAGE : perl perl_blat.pl -BLAT_PORT <BLAT_PORT> -MIN_SCORE <MIN_SCORE> -MIN_IDENTITY <MIN_IDENTITY> -BLAT_PATH <PATH TO BLAT FOLDER> -REF_FILE <PATH TO 2bit file> -INPUT_FILE <INPUT CONFIG FILE> -OUTPUT_FILE <OUTPUT FILE>\n";
}

#parsing the arguments

#unless(-d $outdir)
#{
#	system("mkdir -p $outdir");
#}
$status=`$blatpath/gfServer status localhost $blatport |wc -l`;
chomp($status);
$count = 0;
while($status < 2 )
{
	if($count > 0)
	{
		$blatport = $blatport+int(rand(1000))+1;
	}
	print "Starting the server\n";
	$sys ="$blatpath/gfServer start -canStop localhost $blatport $reffile &";
	print "$sys\n";
	system($sys);
	sleep(300);
	$status=`$blatpath/gfServer status localhost $blatport |wc -l`;
	chomp($status);
	$count++;
	if($count > 5)
	{
		die "something wrong with gfServer or command . Failed 5 times\n";
	}
}	
print "querying \n";
$sys = "$blatpath/gfClient localhost $blatport / $inputfile $outputfile -minScore=$minScore -minIdentity=$minidentity";
print "$sys\n";
system($sys);
print "stoping the server\n";
#$sys = "$blatpath/gfServer stop localhost $blatport";
$pid = `ps|grep gfServer|head -1|cut -f1 -d ' '`;
$sys ="kill -9 $pid";
print "$sys\n";
system($sys);
