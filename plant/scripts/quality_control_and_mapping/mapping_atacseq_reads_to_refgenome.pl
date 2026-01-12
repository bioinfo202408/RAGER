#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;
use File::Find;
use File::Basename;
use File::Spec;

my ($inputdir, $outputdir, $indexdir, $fastqdir, $bowtie2dir, $picarddir, $threads, $help);

GetOptions(
    "inputdir|i=s"   => \$inputdir,
    "outputdir|o=s"  => \$outputdir,
    "indexdir=s"     => \$indexdir,
    "fastqdir=s"     => \$fastqdir,
    "bowtie2dir=s"   => \$bowtie2dir,
    "picarddir=s"    => \$picarddir,
    "threads=s"      => \$threads,
    "help!"          => \$help,
);

# 智能检测样本类型 (单端/双端)
my %samples;
find({
    wanted => sub {
        return unless -f $_;
        return unless /\.fq$/;

        my $full_path = File::Spec->rel2abs($_);
        my $basename = basename($full_path);
        my $dirname  = dirname($full_path);

        # 检测双端数据
        if ($basename =~ /^(.*?)_1_val_1\.fq$/) {
            my $sample_id = $1;
            my $r2_file = File::Spec->catfile($dirname, "${sample_id}_2_val_2.fq");
            if (-e $r2_file) {
                $samples{$sample_id} = {
                    type => 'paired',
                    r1   => $full_path,
                    r2   => $r2_file
                };
            }
        }
        # 检测单端数据
        elsif ($basename =~ /^(.*?)_1_trimmed\.fq$/) {
            my $sample_id = $1;
            $samples{$sample_id} = {
                type => 'single',
                r1   => $full_path
            };
        }
    },
    no_chdir => 1
}, $inputdir);

# 处理每个样本
foreach my $sample_id (keys %samples) {
    my $data_type = $samples{$sample_id}{type};
    my $r1_path = $samples{$sample_id}{r1};
    my $r2_path = $samples{$sample_id}{r2} if $data_type eq 'paired';

    # 创建输出目录
    unless (-d "$outputdir/$bowtie2dir/$sample_id") {
        mkpath("$outputdir/$bowtie2dir/$sample_id", 0644);
        if ($@) {
            die "创建目录失败: $outputdir/$bowtie2dir/$sample_id\n$@";
        }
    }

    # 生成处理脚本
    open(my $sh_fh, ">", "$outputdir/$bowtie2dir/$sample_id/${sample_id}_expcal.sh") 
        or die "无法创建脚本文件: $!";

    # Bowtie2比对命令
    unless (-e "$outputdir/$bowtie2dir/$sample_id/$sample_id.bam") {
        print $sh_fh "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 \\\n";
        print $sh_fh "        --no-mixed --no-discordant \\\n";
        print $sh_fh "        --rg-id $sample_id --rg SM:$sample_id \\\n";

        if ($data_type eq 'paired') {
            print $sh_fh "        -1 $r1_path -2 $r2_path \\\n";
        } else {
            print $sh_fh "        -U $r1_path \\\n";
        }

        print $sh_fh "        -S $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam 2> $outputdir/$bowtie2dir/$sample_id/mapping_summary.txt\n\n";
    }

    # SAM处理流程
    unless (-e "$outputdir/$bowtie2dir/$sample_id/$sample_id.bam") {
        print $sh_fh "samtools view -H $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam > ";
        print $sh_fh "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam\n";
        print $sh_fh "grep 'AS:' $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam | ";
        print $sh_fh "grep -v 'XS:' >> $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam\n";
        print $sh_fh "samtools view -bS $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam \\\n";
        print $sh_fh "        -o $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.bam\n\n";
    }

    # 排序和去重
    unless (-e "$outputdir/$bowtie2dir/$sample_id/$sample_id.bam") {
        print $sh_fh "samtools sort -@ $threads \\\n";
        print $sh_fh "        -o $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.bam \\\n";
        print $sh_fh "        $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.bam\n\n";

        print $sh_fh "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates \\\n";
        print $sh_fh "        I=$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.bam \\\n";
        print $sh_fh "        O=$outputdir/$bowtie2dir/$sample_id/$sample_id.bam \\\n";
        print $sh_fh "        METRICS_FILE=$outputdir/$bowtie2dir/$sample_id/${sample_id}.metricsFile \\\n";
        print $sh_fh "        VALIDATION_STRINGENCY=LENIENT \\\n";
        print $sh_fh "        REMOVE_DUPLICATES=true \\\n";
        print $sh_fh "        ASSUME_SORT_ORDER=coordinate\n\n";
        
        print $sh_fh "samtools index $outputdir/$bowtie2dir/$sample_id/$sample_id.bam\n";
        print $sh_fh "rm -f $outputdir/$bowtie2dir/$sample_id/accepted_hits{.sam,_NHi1.sam,_NHi1.bam,_NHi1.sorted.bam}\n\n";
    }

    close $sh_fh;
}
