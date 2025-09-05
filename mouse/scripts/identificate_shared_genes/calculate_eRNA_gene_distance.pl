#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($enhfile, $genefile, $outfile, $help);
GetOptions(
    "enhfile=s" => \$enhfile,
    "genefile=s" => \$genefile,
    "outfile=s" => \$outfile,
    "help!" => \$help,
);

sub min {
    my $mn = $_[0];
    for my $e (@_) { 
        $mn = $e if ($e < $mn);
    }
    return $mn;
}

my %enhHash;
open(ENH, "<$enhfile") or die "$!\n";
while (<ENH>) {
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\t/, $line;
    $enhHash{$fieldValues[3]}{"chrom"} = $fieldValues[0];
    $enhHash{$fieldValues[3]}{"start"} = $fieldValues[1];
    $enhHash{$fieldValues[3]}{"end"} = $fieldValues[2];
}
close ENH;

my %geneposHash;
open(GENE, "<$genefile") or die "$!\n";
while (<GENE>) {
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\t/, $line;
    my $geneid = $fieldValues[3];
    my $genestart = $fieldValues[1];
    my $geneend = $fieldValues[2];

    $geneposHash{$fieldValues[0]}{$geneid}{"start"} = $genestart;
    $geneposHash{$fieldValues[0]}{$geneid}{"end"} = $geneend;
}
close GENE;

open(OUT, ">$outfile") or die "$!\n";
foreach my $enhid (keys %enhHash) {
    my $enhchrom = $enhHash{$enhid}{"chrom"};
    my $enhstart = $enhHash{$enhid}{"start"};
    my $enhend = $enhHash{$enhid}{"end"};
    my @geneids = keys %{$geneposHash{$enhchrom}};
    foreach my $geneid (@geneids) {
        my $genestart = $geneposHash{$enhchrom}{$geneid}{"start"};
        my $geneend = $geneposHash{$enhchrom}{$geneid}{"end"};
        my $min = min(abs($geneend - $enhstart), abs($genestart - $enhend));
        if ($min < 1000000) {
            print OUT "$enhid\t$geneid\t$min\n";
        }
    }
}
close OUT;
