#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($motifdir, $outfile, $help);

GetOptions(
    "motifdir=s" => \$motifdir,
    "outfile=s"  => \$outfile,
    "help!"      => \$help,
);

my @motifseqfiles = `find $motifdir -name "sequences.tsv"`;

open(my $out_fh, '>', $outfile) or die "Cannot open $outfile: $!";

foreach my $motifseqfile (@motifseqfiles) {
    chomp $motifseqfile;

    (my $enfile = $motifseqfile) =~ s/sequences/ame/;
    next unless -e $enfile;

    open(my $en_fh, '<', $enfile) or die "Cannot open $enfile: $!";
    my ($stat, $adj_pvalue, $motif_Name) = ('No', '', '');
    while (<$en_fh>) {
        chomp;
        next if /^rank/;
        my @fieldValues = split /\t/;
        if (defined $fieldValues[6] && $fieldValues[6] < 0.05) {
            $stat = 'Yes';
        }
        $adj_pvalue = $fieldValues[6] if defined $fieldValues[6];
        $motif_Name  = $fieldValues[3] if defined $fieldValues[3];
        last;
    }
    close $en_fh;

open(my $ms_fh, '<', $motifseqfile) or die "Cannot open $motifseqfile: $!";
my %seen;

while (<$ms_fh>) {
    chomp;
    my @fieldValues = split /\s+/;
    if (defined $fieldValues[6]
        && $fieldValues[6] eq 'tp'
        && $fieldValues[3] !~ /shuf/) {
        
        next if exists $seen{$fieldValues[3]};

        $seen{$fieldValues[3]} = 1;

        print $out_fh join("\t",
            $fieldValues[1],
            $motif_Name,
            $fieldValues[3],
            $fieldValues[5],
            $stat,
            $adj_pvalue
        ), "\n";
    }
}
close $ms_fh;

}

close $out_fh;
