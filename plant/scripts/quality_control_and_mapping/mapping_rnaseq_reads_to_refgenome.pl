#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use warnings;
use File::Path qw(mkpath);
use File::Find;
use File::Basename;
use File::Spec;

my ($inputdir,$outputdir,$indexdir,$hisat2dir,$stringtiedir,$threads,$help,$gfffile,$picarddir);

GetOptions(
    "inputdir|i=s" => \$inputdir,
    "outputdir|o=s" => \$outputdir,
    "indexdir=s" => \$indexdir,
    "hisat2dir=s" => \$hisat2dir,
    "stringtiedir=s" => \$stringtiedir,
    "picarddir=s" => \$picarddir,
    "threads=s" => \$threads,
    "help!" => \$help,
    "gfffile|g=s" => \$gfffile,
);

# 自动检测样本类型（单端/双端）
my %samples;
find({ 
    wanted => sub {
        return unless -f $_;
        return unless /\.fq$/;
        
        my $full_path = $File::Find::name;
        my $basename = basename($full_path);
        my $dirname = dirname($full_path);
        
        # 检测双端数据
        if ($basename =~ /^(.*?)_1_val_1\.fq$/) {
            my $sample_id = $1;
            my $r2 = File::Spec->catfile($dirname, "${sample_id}_2_val_2.fq");
            if (-e $r2) {
                $samples{$sample_id} = {
                    type => 'paired',
                    r1 => $full_path,
                    r2 => $r2
                };
            }
        }
        # 检测单端数据
        elsif ($basename =~ /^(.*?)_1_trimmed\.fq$/) {
            my $sample_id = $1;
            $samples{$sample_id} = {
                type => 'single',
                r1 => $full_path
            };
        }
    },
    no_chdir => 1
}, $inputdir);

# 处理每个样本
foreach my $sample_id (keys %samples) {
    my $type = $samples{$sample_id}{type};
    my $r1 = $samples{$sample_id}{r1};
    my $r2 = $samples{$sample_id}{r2} if $type eq 'paired';

    # 创建输出目录
    mkpath("$outputdir/$hisat2dir/$sample_id", 0644) unless -d "$outputdir/$hisat2dir/$sample_id";
    mkpath("$outputdir/$stringtiedir/$sample_id", 0644) unless -d "$outputdir/$stringtiedir/$sample_id";

    # 生成处理脚本
    open(SH, ">$outputdir/$hisat2dir/$sample_id/${sample_id}_expcal.sh") or die "$!\n";

    # HISAT2比对命令
    if (!-e "$outputdir/$hisat2dir/$sample_id/$sample_id.bam") {
        if ($type eq 'paired') {
            print SH "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id \\\n";
            print SH "      -1 $r1 -2 $r2 \\\n";
        } else {
            print SH "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id \\\n";
            print SH "      -U $r1 \\\n";
        }
        print SH "      --summary-file $outputdir/$hisat2dir/$sample_id/mapping_summary.txt \\\n";
        print SH "      -S $outputdir/$hisat2dir/$sample_id/accepted_hits.sam\n\n";
    }

    # 后续处理步骤（保持原有逻辑不变）
    if (!-e "$outputdir/$hisat2dir/$sample_id/$sample_id.bam") {
        print SH "grep -v -E -w 'NH:i:[2-9]|NH:i:1[0-9]' $outputdir/$hisat2dir/$sample_id/accepted_hits.sam > $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sam\n";
        print SH "samtools view -bS $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sam -o $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
        print SH "samtools sort -@ $threads -o $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
        print SH "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates \\\n";
        print SH "      I=$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam \\\n";
        print SH "      O=$outputdir/$hisat2dir/$sample_id/$sample_id.bam \\\n";
        print SH "      METRICS_FILE=$outputdir/$hisat2dir/$sample_id/${sample_id}.metricsFile \\\n";
        print SH "      VALIDATION_STRINGENCY=LENIENT \\\n";
        print SH "      REMOVE_DUPLICATES=true \\\n";
        print SH "      ASSUME_SORT_ORDER=coordinate\n\n";
        print SH "samtools index $outputdir/$hisat2dir/$sample_id/$sample_id.bam\n";
        print SH "rm -f $outputdir/$hisat2dir/$sample_id/accepted_hits{.sam,_NHi1.sam,_NHi1.bam,_NHi1.sorted.bam}\n\n";
    }

    # StringTie定量分析
    if (!-e "$outputdir/$stringtiedir/$sample_id/transcripts.gtf") {
        print SH "stringtie -p $threads -e -B \\\n";
        print SH "          -G $gfffile \\\n";
        print SH "          -A $outputdir/$stringtiedir/$sample_id/gene_abund.tab \\\n";
        print SH "          -o $outputdir/$stringtiedir/$sample_id/transcripts.gtf \\\n";
        print SH "          $outputdir/$hisat2dir/$sample_id/$sample_id.bam\n";
    }

    close SH;
}
