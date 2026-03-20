#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($fastafile,$outputdir,$motifdir,$help);
GetOptions(
	"fastafile=s" => \$fastafile,
    "outputdir=s" => \$outputdir,
	"motifdir=s" => \$motifdir,
	"help!" => \$help,
);

my @motiffiles = `find $motifdir -name "*meme"`;
foreach my $motiffile (@motiffiles){
	$motiffile =~ /.*\/(.*)\.meme/;
	my $motifname = $1;
	if(!-e "$outputdir/$motifname"){
		mkpath("$outputdir/$motifname",0644);
		if($@){
			print "Make path $outputdir/$motifname failed:$@\n";
			exit(1);
		}

		my $out = system("ame --control --shuffle-- --oc $outputdir/$motifname $fastafile $motiffile");
		if($out == 0){
			print "The task of $motifname is successfully submitted\n";
		}
	}
}
