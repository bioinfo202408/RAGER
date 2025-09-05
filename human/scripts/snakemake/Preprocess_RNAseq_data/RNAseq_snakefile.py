from snakemake.io import directory, ancient

# Load configuration file
configfile: "config.yaml"

# Extract sample information from config
SAMPLES = {}
for sample, info in config["samples"].items():
    SAMPLES[sample] = info["type"]

# Load paths from config
BASE_DIR = config["paths"]["base_dir"]
SCRIPTS_DIR = config["paths"]["scripts_dir"]
RNAseq_BASE_DIR = config["paths"]["rnaseq_base_dir"]
RNAseq_FASTQ_DIR = config["paths"]["fastq_dir"]
RNAseq_QC_DIR = config["paths"]["qc_dir"]
RNAseq_HISAT2_DIR = config["paths"]["hisat2_dir"]
RNAseq_STRINGTIE_DIR = config["paths"]["stringtie_dir"]
RNAseq_INDEX_DIR = config["paths"]["index_dir"]

# Load reference files from config
GFFFILE = config["references"]["gtf_file"]
RseQC_BEDFILE = config["references"]["bed_file"]

# Generate final files list based on sample types
FINAL_FILES = []
for sample, stype in SAMPLES.items():
    if stype == "paired":
        FINAL_FILES.extend([
            f"{RNAseq_QC_DIR}/{sample}_1_val_1.fq",
            f"{RNAseq_QC_DIR}/{sample}_2_val_2.fq"
        ])
    else:
        FINAL_FILES.append(f"{RNAseq_QC_DIR}/{sample}_1_trimmed.fq")

rule all:
    input:
        FINAL_FILES,
        f"{RNAseq_QC_DIR}/multiqc_report.html",
        expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}_expcal.sh", sample=SAMPLES),
        expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}.bam", sample=SAMPLES),
        expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}.bam.bai", sample=SAMPLES),
        expand(f"{RNAseq_STRINGTIE_DIR}/{{sample}}/gene_abund.tab", sample=SAMPLES),
        expand(f"{RNAseq_STRINGTIE_DIR}/{{sample}}/transcripts.gtf", sample=SAMPLES),
        expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/mapping_summary.txt", sample=SAMPLES),
        f"{RNAseq_BASE_DIR}/hisat2file.geneBodyCoverage.curves.pdf",
        f"{RNAseq_BASE_DIR}/hisat2file.geneBodyCoverage.heatMap.pdf",
        f"{RNAseq_HISAT2_DIR}/RNAseq_alignment_rates_plot.pdf",
        f"{RNAseq_STRINGTIE_DIR}/gene_TPM.txt",
        f"{RNAseq_STRINGTIE_DIR}/RNAseq_sample_clustering.pdf",
        f"{RNAseq_STRINGTIE_DIR}/3DPCA.pdf",
        f"{RNAseq_STRINGTIE_DIR}/samplelist.txt",
        f"{RNAseq_STRINGTIE_DIR}/gene_count_matrix.count",
        f"{RNAseq_STRINGTIE_DIR}/{config['params']['deseq2']['treatment_group']}_vs_{config['params']['deseq2']['control_group']}_DEG.csv"

rule paired_trimmed:
    input:
        p1 = f"{RNAseq_FASTQ_DIR}/{{sample}}_1.fastq",
        p2 = f"{RNAseq_FASTQ_DIR}/{{sample}}_2.fastq"
    output:
        qc1 = f"{RNAseq_QC_DIR}/{{sample}}_1_val_1.fq",
        qc2 = f"{RNAseq_QC_DIR}/{{sample}}_2_val_2.fq"    
    params:
        output_dir = RNAseq_QC_DIR,
        phred = config["params"]["trim_galore"]["phred"],
        quality = config["params"]["trim_galore"]["quality"],
        length = config["params"]["trim_galore"]["length_paired"],
        stringency = config["params"]["trim_galore"]["stringency_paired"]
    shell:
        """
        trim_galore --phred{params.phred} -q {params.quality} --length {params.length} \
        --stringency {params.stringency} --paired --fastqc -o {params.output_dir} {input.p1} {input.p2}
        """

rule single_trimmed:
    input:
        p1 = f"{RNAseq_FASTQ_DIR}/{{sample}}_1.fastq"
    output:
        qc1 = f"{RNAseq_QC_DIR}/{{sample}}_1_trimmed.fq"  
    params:
        output_dir = RNAseq_QC_DIR,
        phred = config["params"]["trim_galore"]["phred"],
        quality = config["params"]["trim_galore"]["quality"],
        length = config["params"]["trim_galore"]["length_single"],
        stringency = config["params"]["trim_galore"]["stringency_single"]
    shell:
        """
        trim_galore --phred{params.phred} -q {params.quality} --length {params.length} \
        --stringency {params.stringency} --fastqc -o {params.output_dir} {input.p1}
        """

rule multiqc:
    input:
        FINAL_FILES
    output:
        f"{RNAseq_QC_DIR}/multiqc_report.html"
    shell:
        """
        multiqc {RNAseq_QC_DIR} -o {RNAseq_QC_DIR}
        """

rule mapping_reads_to_genome:
    input:
        lambda wildcards: (
            [
                f"{RNAseq_QC_DIR}/{wildcards.sample}_1_val_1.fq",
                f"{RNAseq_QC_DIR}/{wildcards.sample}_2_val_2.fq"
            ] if SAMPLES[wildcards.sample] == "paired" else
            [
                f"{RNAseq_QC_DIR}/{wildcards.sample}_1_trimmed.fq"
            ]
        )
    output:
        result = f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}_expcal.sh"
    params:
        hisat2_dir = config["output_dirs"]["hisat2_subdir"],
        stringtie_dir = config["output_dirs"]["stringtie_subdir"],
        threads = config["params"]["hisat2"]["threads"],
        index_name = config["params"]["hisat2"]["index_name"]
    shell:
        f"""
        perl {SCRIPTS_DIR}/quality_control_and_mapping/mapping_rnaseq_reads_to_refgenome.pl \
            --inputdir {input} \
            --outputdir {RNAseq_BASE_DIR} \
            --indexdir {RNAseq_INDEX_DIR}/{params.index_name} \
            --picarddir {BASE_DIR} \
            --hisat2dir {params.hisat2_dir} \
            --stringtiedir {params.stringtie_dir} \
            --gfffile {GFFFILE} \
            --threads {params.threads}
        """

# Run the generated shell script
rule run_expcal_sh:
    input:
        sh = f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}_expcal.sh"
    output:
        bam = f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}.bam",
        bai = f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}.bam.bai",
        TPM = f"{RNAseq_STRINGTIE_DIR}/{{sample}}/gene_abund.tab",
        gtf = f"{RNAseq_STRINGTIE_DIR}/{{sample}}/transcripts.gtf",
        summaries = f"{RNAseq_HISAT2_DIR}/{{sample}}/mapping_summary.txt"
    shell:
        "bash {input.sh}"

rule gene_body_coverage:
    input:
        bed = config["references"]["bed_file"],
        bams = expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/{{sample}}.bam", sample=SAMPLES)
    output:
        pdf1 = f"{RNAseq_BASE_DIR}/hisat2file.geneBodyCoverage.curves.pdf",
        pdf2 = f"{RNAseq_BASE_DIR}/hisat2file.geneBodyCoverage.heatMap.pdf"
    params:
        output_prefix = f"{RNAseq_BASE_DIR}/hisat2file"
    run:
        bam_files = ",".join(input.bams)
        shell_cmd = f"geneBody_coverage.py -r {input.bed} -i {bam_files} -o {params.output_prefix}"
        shell(shell_cmd)

rule plot_alignment_rates:
    input:
        summaries = expand(f"{RNAseq_HISAT2_DIR}/{{sample}}/mapping_summary.txt", sample=SAMPLES)
    output:
        plot = f"{RNAseq_HISAT2_DIR}/RNAseq_alignment_rates_plot.pdf"
    params:
        script_dir = config["paths"]["scripts_dir"]
    shell:
        """
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/RNAseq_alignment_rates.R "{input.summaries}" {output.plot}
        """

rule extract_expression_data:
    input:
        gene_abund = expand(f"{RNAseq_STRINGTIE_DIR}/{{sample}}/gene_abund.tab", sample=SAMPLES)
    output:
        tpm = f"{RNAseq_STRINGTIE_DIR}/gene_TPM.txt"
    params:
        script_dir = config["paths"]["scripts_dir"],
        stringtie_dir = RNAseq_STRINGTIE_DIR
    shell:
        """
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/extr_expr_data.R {params.stringtie_dir} {output.tpm}
        """

rule generate_cluster_plots:
    input:
        TPM = f"{RNAseq_STRINGTIE_DIR}/gene_TPM.txt"
    output:
        p = f"{RNAseq_STRINGTIE_DIR}/RNAseq_sample_clustering.pdf"
    params:
        script_dir = config["paths"]["scripts_dir"],
        # Extract sample names from config
        sample_ids = ",".join([f'"{sample}"' for sample in config["samples"].keys()]),
        # Extract sample labels from config  
        sample_labels = ",".join([f'"{info["label"]}"' for info in config["samples"].values()])
    shell:
        """
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/RNAseq_cluster_analysis.R \
            {input.TPM} \
            {params.sample_ids} \
            {params.sample_labels} \
            {output.p}
        """

rule generate_pca_plots:
    input:
        TPM = f"{RNAseq_STRINGTIE_DIR}/gene_TPM.txt"
    output:
        pc1 = f"{RNAseq_STRINGTIE_DIR}/2DPCA_PC1_PC2.pdf",
        pc2 = f"{RNAseq_STRINGTIE_DIR}/2DPCA_PC1_PC3.pdf", 
        pc3 = f"{RNAseq_STRINGTIE_DIR}/2DPCA_PC2_PC3.pdf",
        pca3d = f"{RNAseq_STRINGTIE_DIR}/3DPCA.pdf"
    params:
        script_dir = config["paths"]["scripts_dir"],
        # Extract sample names from config
        sample_names = ",".join([f'"{sample}"' for sample in config["samples"].keys()]),
        # Extract sample labels from config
        sample_labels = ",".join([f'"{info["label"]}"' for info in config["samples"].values()]),
        # Extract group names from config
        group_names = ",".join([f'"{group}"' for group in config["params"]["analysis_groups"]["group_names"]]),
        # Extract sample replicates from config
        sample_replicates = ",".join([str(rep) for rep in config["params"]["analysis_groups"]["sample_replicates"]])
    shell:
        """
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/RNAseq_PCA_analysis.R \
            {input.TPM} \
            {params.sample_names} \
            {params.sample_labels} \
            {params.group_names} \
            {params.sample_replicates} \
            {output.pc1} {output.pc2} {output.pc3} {output.pca3d}
        """

rule construct_samplelist:
    input:
         gtf= expand(f"{RNAseq_STRINGTIE_DIR}/{{sample}}/transcripts.gtf", sample=SAMPLES)
    output:
        list = f"{RNAseq_STRINGTIE_DIR}/samplelist.txt"
    params:
        script_dir = config["paths"]["scripts_dir"],
        # Extract sample labels from config in order
        sample_labels = ",".join([f'"{config["samples"][sample]["label"]}"' for sample in SAMPLES.keys()])
    run:
        gtf_files = ",".join(input.gtf)
        shell_cmd = f"""
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/construct_samplelist.R \
            {params.sample_labels} \
            {gtf_files} \
            {output.list}
        """
        shell(shell_cmd)

rule prepDE:
    input:
        list = f"{RNAseq_STRINGTIE_DIR}/samplelist.txt"
    output:
        count1 = f"{RNAseq_STRINGTIE_DIR}/gene_count_matrix.count",
        count2 = f"{RNAseq_STRINGTIE_DIR}/transcript_count_matrix.count"
    params:
        script_dir = config["paths"]["scripts_dir"]
    shell:
        """
        python3 {params.script_dir}/quality_control_and_mapping/prepDE.py3 -i {input.list} -g {output.count1} -t {output.count2}
        """

rule DEseq2:
    input:
        count_file = f"{RNAseq_STRINGTIE_DIR}/gene_count_matrix.count"
    output:
        DEexpr = f"{RNAseq_STRINGTIE_DIR}/{config['params']['deseq2']['treatment_group']}_vs_{config['params']['deseq2']['control_group']}_DEG.csv"
    params:
        script_dir = config["paths"]["scripts_dir"],
        # Extract group information for each sample from config
        sample_groups = ",".join([f'"{config["samples"][sample]["group"]}"' for sample in SAMPLES.keys()]),
        # Extract control and treatment groups from config
        control_group = config["params"]["deseq2"]["control_group"],
        treatment_group = config["params"]["deseq2"]["treatment_group"],
        # Create comparison groups parameter
        comparison_groups = f'"{config["params"]["deseq2"]["control_group"]}","{config["params"]["deseq2"]["treatment_group"]}"'
    shell:
        """
        Rscript --vanilla {params.script_dir}/quality_control_and_mapping/DEseq2.R \
            {input.count_file} \
            {params.sample_groups} \
            {params.comparison_groups} \
            {output.DEexpr}
        """