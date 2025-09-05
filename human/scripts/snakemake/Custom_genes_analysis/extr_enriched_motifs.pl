#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($tfmapfile, $motifdir, $outfile, $help);
GetOptions(
    "tfmapfile=s" => \$tfmapfile,
    "motifdir=s" => \$motifdir,
    "outfile=s" => \$outfile,
    "help!" => \$help,
);

my %tfmapHash;
open(MAP, "<$tfmapfile") or die "$!\n";
while (<MAP>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^Homo_sapiens/) {
        my @fieldValues = split /\t/, $line;
        $tfmapHash{uc($fieldValues[1])} = $fieldValues[2];
    }
}
close MAP;

my @enfiles = `find $motifdir -name "ame.tsv"`;
open(OUT, ">$outfile") or die "$!\n";
foreach my $enfile (@enfiles) {
    chomp $enfile;
    open(EN, "<$enfile") or die "$!\n";
    while (<EN>) {
        my $line = $_;
        chomp $line;
        next if ($line =~ /^rank/);
        my @fieldValues = split /\t/, $line;

        if (defined $fieldValues[6] && $fieldValues[6] < 0.05) {
            if ($fieldValues[3] =~ /\::/) {
                my @parts = split /::/, $fieldValues[3];
                foreach my $part (@parts) {
                    if (exists $tfmapHash{uc($part)}) {
                        my $tfid = $tfmapHash{uc($part)};
                        my $tfname = uc($part);
                        print OUT "$tfid\t$tfname\t$fieldValues[6]\n";
                    }
                }
            } else {
                if (exists $tfmapHash{uc($fieldValues[3])}) {
                    my $tfid = $tfmapHash{uc($fieldValues[3])};
                    my $tfname = uc($fieldValues[3]);
                    print OUT "$tfid\t$tfname\t$fieldValues[6]\n";
                }
            }
        }
    }
    close EN;
}
close OUT;
#perl extr_enriched_motifs.pl --tfmapfile /home/yjliu/RAGER/geneanno/Mus_musculus_TF.txt --motifdir /home/yjliu/RAGER/dataset/mouse_embryonic_sc/shared_UPgene_motifenrichdata/ --outfile /home/yjliu/RAGER/dataset/mouse_embryonic_sc/shared_UPgene_motifenrichdata/enriched_tfs_list.txt