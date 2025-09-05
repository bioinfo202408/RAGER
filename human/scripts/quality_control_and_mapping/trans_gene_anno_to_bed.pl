#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($infile,$outfile,$help);

GetOptions(
	"infile|i=s" => \$infile,
	"outfile|o=s" => \$outfile,
	"help!" => \$help,
);

open IN, "<$infile" or die "$!\n";
open(OUT,">$outfile") or die "$!\n";
while(<IN>){
	my $line = $_;
	chomp $line;
	next if($line =~ /^#/);
	my @fieldValues = split /\t/,$line;
	if($fieldValues[2] eq "gene"){
		$fieldValues[8] =~ /gene_id\s+"([^"]+)"/;
		my $geneid = $1;
		print OUT "$fieldValues[0]\t$fieldValues[3]\t$fieldValues[4]\t$geneid\n";
	}
}
close IN;
close OUT;
#perl /home/yjliu/RAGER/case2_test/scripts/qc_mapping/trans_gene_anno_to_bed.pl --infile /home/yjliu/RAGER/case2_test/reference_annotation_file/geneanno/gencode.v44.annotation.gtf --outfile /home/yjliu/RAGER/case2_test/reference_annotation_file/geneanno/gencode.v44.annotation.bed