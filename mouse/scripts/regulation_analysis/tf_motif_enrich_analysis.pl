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

# perl tf_motif_enrich_analysis.pl --fastafile /data_p2/cjxProj/geneanno/promoter_seqs.fa --outputdir /data_p2/cjxProj/datasets/motifenrichdata --motifdir /data_p2/cjxProj/geneanno/JASPAR2024_CORE

