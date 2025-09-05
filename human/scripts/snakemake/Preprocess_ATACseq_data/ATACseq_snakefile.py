# Load configuration
configfile: "config.yaml"

# Extract samples information from config
SAMPLES = {sample: info["type"] for sample, info in config["samples"].items()}

# Extract paths from config
BASE_DIR = config["paths"]["base_dir"]
SCRIPTS_DIR = config["paths"]["scripts_dir"]
ATACseq_BASE_DIR = config["paths"]["atacseq_base_dir"]
ATACseq_FASTQ_DIR = config["paths"]["fastq_dir"]
ATACseq_QC_DIR = config["paths"]["qc_dir"]
ATACseq_BOWTIE2_DIR = config["paths"]["bowtie2_dir"]
ATACseq_INDEX_DIR = config["paths"]["index_dir"]
ATACseq_MACS2_DIR = config["paths"]["macs2_dir"]
GENEANNO_DIR = config["paths"]["geneanno_dir"]

# Generate final files list
FINAL_FILES = []
for sample, stype in SAMPLES.items():
    if stype == "paired":
        FINAL_FILES.extend([
            f"{ATACseq_QC_DIR}/{sample}_1_val_1.fq",
            f"{ATACseq_QC_DIR}/{sample}_2_val_2.fq"
        ])
    else:
        FINAL_FILES.append(f"{ATACseq_QC_DIR}/{sample}_1_trimmed.fq")

rule all:
    input:
        FINAL_FILES,
        f"{ATACseq_QC_DIR}/multiqc_report.html",
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}_expcal.sh", sample=SAMPLES),
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam", sample=SAMPLES),
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam.bai", sample=SAMPLES),
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/mapping_summary.txt", sample=SAMPLES),
        f"{ATACseq_BOWTIE2_DIR}/ATACseq_alignment_rates_plot.pdf",
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/fragmentSizeDistribution.pdf", sample=SAMPLES),
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/heatmap.pdf", sample=SAMPLES),
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/coverage_curve_plot.pdf", sample=SAMPLES),
        f"{config['references']['bed_file']}",
        expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}.bw", sample=SAMPLES),
        f"{ATACseq_BOWTIE2_DIR}/all_scaled.tab",
        f"{ATACseq_BOWTIE2_DIR}/sample_clustering_reorder.pdf",
        f"{ATACseq_BOWTIE2_DIR}/3DPCA.pdf",
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam",
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam.bai",
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam",
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam.bai",
        f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_summits.bed",
        f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_summits.bed",
        f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_peak.csv",
        f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_peak.csv",
        f"{ATACseq_MACS2_DIR}/peaks_diverging_bar.pdf",
        f"{ATACseq_MACS2_DIR}/all_peak_bar.pdf"

rule paired_trimmed:
    input:
        p1 = f"{ATACseq_FASTQ_DIR}/{{sample}}_1.fastq",
        p2 = f"{ATACseq_FASTQ_DIR}/{{sample}}_2.fastq"
    output:
        qc1 = f"{ATACseq_QC_DIR}/{{sample}}_1_val_1.fq",
        qc2 = f"{ATACseq_QC_DIR}/{{sample}}_2_val_2.fq"    
    params:
        output_dir=ATACseq_QC_DIR,
        phred=config["params"]["trim_galore"]["phred"],
        quality=config["params"]["trim_galore"]["quality"],
        length=config["params"]["trim_galore"]["length_paired"],
        stringency=config["params"]["trim_galore"]["stringency_paired"]
    shell:
        """
        trim_galore --phred{params.phred} -q {params.quality} --length {params.length} --stringency {params.stringency} --paired --fastqc -o {params.output_dir} {input.p1} {input.p2}
        """

rule single_trimmed:
    input:
        p1 = f"{ATACseq_FASTQ_DIR}/{{sample}}_1.fastq"
    output:
        qc1 = f"{ATACseq_QC_DIR}/{{sample}}_1_trimmed.fq"  
    params:
        output_dir=ATACseq_QC_DIR,
        phred=config["params"]["trim_galore"]["phred"],
        quality=config["params"]["trim_galore"]["quality"],
        length=config["params"]["trim_galore"]["length_single"],
        stringency=config["params"]["trim_galore"]["stringency_single"]
    shell:
        """
        trim_galore --phred{params.phred} -q {params.quality} --length {params.length} --stringency {params.stringency} --fastqc -o {params.output_dir} {input.p1}
        """

rule multiqc:
    input:
        FINAL_FILES
    output:
        f"{ATACseq_QC_DIR}/multiqc_report.html"
    shell:
        """
        multiqc {ATACseq_QC_DIR} -o {ATACseq_QC_DIR}
        """

rule mapping_reads_to_genome:
    input:
        lambda wildcards: (
            [
                f"{ATACseq_QC_DIR}/{wildcards.sample}_1_val_1.fq",
                f"{ATACseq_QC_DIR}/{wildcards.sample}_2_val_2.fq"
            ] if SAMPLES[wildcards.sample] == "paired" else
            [
                f"{ATACseq_QC_DIR}/{wildcards.sample}_1_trimmed.fq"
            ]
        )
    output:
        result=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}_expcal.sh"
    params:
        bowtie2_dir=config["output_dirs"]["bowtie2_subdir"],
        threads=config["params"]["bowtie2"]["threads"],
        index_name=config["references"]["bowtie2_index"]
    shell:
        f"""
        perl {SCRIPTS_DIR}/quality_control_and_mapping/mapping_atacseq_reads_to_refgenome.pl \
            --inputdir {input} \
            --outputdir {ATACseq_BASE_DIR} \
            --indexdir {ATACseq_INDEX_DIR}/{{params.index_name}} \
            --picarddir {BASE_DIR} \
            --bowtie2dir {{params.bowtie2_dir}} \
            --threads {{params.threads}}
        """

rule run_expcal_sh:
    input:
        sh=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}_expcal.sh"
    output:
        bam=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam",
        bai = f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam.bai",
        summaries=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/mapping_summary.txt"
    shell:
        """
        bash {input.sh}
        touch {output.bam} {output.summaries}
        """

rule plot_alignment_rates:
    input:
        summaries=expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}/mapping_summary.txt", sample=SAMPLES)
    output:
        plot=f"{ATACseq_BOWTIE2_DIR}/ATACseq_alignment_rates_plot.pdf"
    shell:
        f"""
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/ATACseq_alignment_rates.R "{input.summaries}" {output.plot}
        """

rule fragmentSize_plots:
    input:
        bam=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam"
    output:
        fragment_size=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/fragmentSizeDistribution.pdf"
    shell:
        f"""
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/fragmentSize_plot.R {input.bam} {output.fragment_size}
        """

rule generate_qc_plots:
    input:
        bam=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam"
    output:
        PT=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/PTscore_plot.pdf",
        NFRscore=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/NFRscore_plot.pdf",
        TSSEscore=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/TSSEscore_plot.pdf",
        cumulative_percentage=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/cumulative_percentage_plot.pdf",
        heatmap=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/heatmap.pdf",
        coverage_curve=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/coverage_curve_plot.pdf"
    shell:
        f"""
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/ATACseqQC.R {input.bam} {ATACseq_BOWTIE2_DIR}/{{sample}}/splited {output.PT} {output.NFRscore} {output.TSSEscore} {output.cumulative_percentage} {output.heatmap} {output.coverage_curve}
        """

rule trans_gene_anno_to_bed:
    input:
        gfffile=config["references"]["gtf_file"]
    output:
        bedfile=config["references"]["bed_file"]
    shell:
        f"""
        perl {SCRIPTS_DIR}/quality_control_and_mapping/trans_gene_anno_to_bed.pl \
            --infile {input.gfffile} \
            --outfile {output.bedfile}
        """

rule bam_to_bw:
    input:
        bam=f"{ATACseq_BOWTIE2_DIR}/{{sample}}/{{sample}}.bam"
    output:
        bw=f"{ATACseq_BOWTIE2_DIR}/{{sample}}.bw"
    params:
        bin_size=config["params"]["bamcoverage"]["bin_size"],
        processors=config["params"]["bamcoverage"]["processors"],
        normalization=config["params"]["bamcoverage"]["normalization"]
    shell:
        """
        bamCoverage -b {input.bam} -of bigwig --binSize {params.bin_size} --ignoreDuplicates --normalizeUsing {params.normalization} --numberOfProcessors {params.processors} -o {output.bw}
        """

rule calATACexpr:
    input:
        bw=expand(f"{ATACseq_BOWTIE2_DIR}/{{sample}}.bw", sample=SAMPLES),
        bedfile=config["references"]["bed_file"]
    output:
        gz=f"{ATACseq_BOWTIE2_DIR}/all_scaled.gz",
        tab=f"{ATACseq_BOWTIE2_DIR}/all_scaled.tab"
    params:
        before=config["params"]["computematrix"]["before"],
        after=config["params"]["computematrix"]["after"],
        bin_size=config["params"]["computematrix"]["bin_size"],
        processors=config["params"]["computematrix"]["processors"]
    run:
        bw_files = " ".join(input.bw)
        shell_cmd = f"computeMatrix reference-point -R {input.bedfile} -S {bw_files} -b {params.before} -a {params.after} --binSize {params.bin_size} -p {params.processors} -o {output.gz} --outFileNameMatrix {output.tab}"
        shell(shell_cmd)

rule sample_clustering_plot:
    input:
        tab = f"{ATACseq_BOWTIE2_DIR}/all_scaled.tab"
    output:
        clustering = f"{ATACseq_BOWTIE2_DIR}/sample_clustering.pdf",
        reordered = f"{ATACseq_BOWTIE2_DIR}/sample_clustering_reorder.pdf"
    params:
        sample_names = ",".join(config["params"]["analysis_groups"]["sample_labels"]),
        col_ranges = ",".join(config["params"]["analysis_groups"]["column_ranges"]),
        order = ",".join(map(str, config["params"]["analysis_groups"]["sample_order"]))
    shell:
        """
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/ATACseq_cluster_analysis.R \
            {input.tab} {output.clustering} {output.reordered} \
            '{params.sample_names}' '{params.col_ranges}' '{params.order}'
        """

rule generate_pca_plots:
    input:
        tab = f"{ATACseq_BOWTIE2_DIR}/all_scaled.tab"
    output:
        pc1   = f"{ATACseq_BOWTIE2_DIR}/2DPCA_PC1_PC2.pdf",
        pc2   = f"{ATACseq_BOWTIE2_DIR}/2DPCA_PC1_PC3.pdf",
        pc3   = f"{ATACseq_BOWTIE2_DIR}/2DPCA_PC2_PC3.pdf",
        pca3d = f"{ATACseq_BOWTIE2_DIR}/3DPCA.pdf"
    params:
        sample_names = ",".join(config["params"]["analysis_groups"]["sample_labels"]),
        col_ranges = ",".join(config["params"]["analysis_groups"]["column_ranges"]),
        group_names = ",".join(config["params"]["analysis_groups"]["group_names"]),
        sample_rep = ",".join(map(str, config["params"]["analysis_groups"]["sample_replicates"])),
    shell:
        """
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/ATACseq_PCA_analysis.R \
            {input.tab} \
            {params.col_ranges}\
            {params.sample_names} \
            {params.group_names} \
            {params.sample_rep} \
            {output.pc1} {output.pc2} {output.pc3} {output.pca3d}
        """

rule merge_treatment_group_bam:
    input:
        bams = [f"{ATACseq_BOWTIE2_DIR}/{sample}/{sample}.bam" for sample in config["params"]["differential_analysis"]["treatment_samples"]]
    output:
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam"
    shell:
        "samtools merge -f {output} {input.bams}"

rule index_treatment_group_bam:
    input:
        bam = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam"
    output:
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam.bai"
    shell:
        "samtools index {input.bam}"

rule merge_control_group_bam:
    input:
        bams = [f"{ATACseq_BOWTIE2_DIR}/{sample}/{sample}.bam" for sample in config["params"]["differential_analysis"]["control_samples"]]
    output:
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam"
    shell:
        "samtools merge -f {output} {input.bams}"

rule index_control_group_bam:
    input:
        bam = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam"
    output:
        f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam.bai"
    shell:
        "samtools index {input.bam}"

rule macs2_callpeak_Up:
    input:
        treatment = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam",
        control   = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam"
    output:
        peakbed = f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_summits.bed"
    params:
        outdir = f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}",
        name   = f"{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}",
        genome = config["params"]["macs2"]["genome_size"],
        q = config["params"]["macs2"]["qvalue"]
    shell:
        """
        macs2 callpeak \
            -t {input.treatment} \
            -c {input.control} \
            -g {params.genome} \
            -n {params.name} \
            --keep-dup all \
            -q {params.q} \
            --outdir {params.outdir}
        """

rule macs2_callpeak_Down:
    input:
        treatment = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['treatment_group']}.bam",
        control   = f"{ATACseq_BOWTIE2_DIR}/{config['params']['differential_analysis']['control_group']}.bam"
    output:
        peakbed = f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_summits.bed"
    params:
        outdir = f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}",
        name   = f"{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}",
        genome = config["params"]["macs2"]["genome_size"],
        q = config["params"]["macs2"]["qvalue"]
    shell:
        """
        macs2 callpeak \
            -t {input.control} \
            -c {input.treatment} \
            -g {params.genome} \
            -n {params.name} \
            --keep-dup all \
            -q {params.q} \
            --outdir {params.outdir}
        """

rule annoATACPeaks:
    input:
        up_peak=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_summits.bed",
        down_peak=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_summits.bed"
    output:
        anno_up=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_peak.csv",
        anno_down=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_peak.csv"
    shell:
        f"""
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/annoATACPeaks.R {input.up_peak} {input.down_peak} {output.anno_up} {output.anno_down}
        """

rule peaks_bar_plots:
    input:
        up_peak=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['treatment_group']}_vs_{config['params']['differential_analysis']['control_group']}_peak.csv",
        down_peak=f"{ATACseq_MACS2_DIR}/{config['params']['differential_analysis']['control_group']}_vs_{config['params']['differential_analysis']['treatment_group']}_peak.csv"
    output:
        plot1=f"{ATACseq_MACS2_DIR}/peaks_diverging_bar.pdf",
        plot2=f"{ATACseq_MACS2_DIR}/all_peak_bar.pdf"
    shell:
        f"""
        Rscript --vanilla {SCRIPTS_DIR}/quality_control_and_mapping/peaks_number_barplot.R {input.up_peak} {input.down_peak} {output.plot1} {output.plot2}
        """
