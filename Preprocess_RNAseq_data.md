# **RAGER table of contents**
1. [Quick start](https://github.com/bioinfo202408/RAGER/blob/main/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_RNAseq_data.md)
3. [Preprocess ATACseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_ATACseq_data.md) 
4. [Joint analysis](https://github.com/bioinfo202408/RAGER/blob/main/Joint_analysis.md)
5. [Custom analysis](https://github.com/bioinfo202408/RAGER/blob/main/Custom_analysis.md)

## **Quick Start**

Welcome to RAGER! This section will guide through the essential steps to configure and run RNA-seq analysis pipeline.

### **Step 1: Open the Configuration File**

Open the `config.yaml` file in the species-specific directory:

```bash
cd ~/PROJECT/RAGER  #Change the directory to PROJECT/RAGER as the current directory
#mouse
vim ./mouse/scripts/snakemake/Preprocess_RNAseq_data/config.yaml  
#human
vim ./human/scripts/snakemake/Preprocess_RNAseq_data/config.yaml 
# you can use any text editor 
```

### **Step 2:Modifying the Configuration File**

The `config.yaml` file contains all parameters needed to customize RAGER for specific analysis. Below is a detailed explanation of each section:

#### **Paths Section**
This section defines the directory structure for analysis:

```yaml
paths:
  base_dir: "datasets"                    # Main data directory
  scripts_dir: "scripts"                  # Directory containing analysis scripts
  rnaseq_base_dir: "datasets/RNAseq"     # RNA-seq analysis base directory
  fastq_dir: "datasets/RNAseq/fastqfile" # Raw FASTQ files location
  qc_dir: "datasets/RNAseq/quality_control_file"  # Quality control outputs
  hisat2_dir: "datasets/RNAseq/hisat2file"        # Alignment outputs
  stringtie_dir: "datasets/RNAseq/stringtiefile"  # Expression quantification outputs
  index_dir: "reference/hisat2index"              # Reference genome index
```

**How to modify**: All directories exist or will be created by the pipeline,it is recommended to use the default path above.

#### **References Section**
Specify the reference files for target organism:

```yaml
#mouse
references:
  gtf_file: "reference/geneanno/gencode.vM25.annotation.gtf"  # Gene annotation file
  bed_file: "reference/mm10_GENCODE_VM25_basic.bed"           # RSeQC software requires gene annotation files 
#human
references:
  gtf_file: "reference/geneanno/gencode.v44.annotation.gtf"  # Gene annotation file
  bed_file: "reference/hg38_GENCODE_V44_Basic.bed"           # RSeQC software requires gene annotation files
```

**How to modify**: 
- Replace with species-specific GTF and BED files

#### **Samples Section**
Define experimental samples and their properties:

```yaml
samples:
  SRR4032350:           # Sample ID (must match FASTQ file prefix)
    type: "paired"       # "paired" or "single" end sequencing
    group: "RNA_A"     # Experimental group/condition
    label: "RNA_A1"    # Human-readable sample label
  SRR4032351: 
    type: "paired"
    group: "RNA_A"
    label: "RNA_A2"
  SRR4032352: 
    type: "paired"
    group: "RNA_N"
    label: "RNA_N1"
  SRR4032353: 
    type: "paired"
    group: "RNA_N"
    label: "RNA_N2"
```

**How to modify**:
- Replace sample IDs with actual sample names
- Set `type` to `"paired"` for paired-end or `"single"` for single-end data
- Assign meaningful `group` names for experimental conditions
- Provide descriptive `label` names for easier interpretation

#### **Parameters Section**
**trim galore parameters:**
```yaml
params:
  trim_galore:                    # Quality control parameters
    phred: 33                     # Quality score encoding
    quality: 25                   # Quality threshold for trimming
    length_paired: 35             # Minimum read length (paired-end)
    length_single: 30             # Minimum read length (single-end)
```
How to modify:
- Adjust quality parameters based on data quality requirements.For details, please refer to [trim galore website](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)


**mapping parameters:**
```yaml
#mouse
  hisat2:                         # Alignment parameters
    threads: 5                    # Number of CPU threads
    index_name: "GRCm38"         # Reference genome index name
#human
  hisat2:                         # Alignment parameters
    threads: 5                    # Number of CPU threads
    index_name: "GRCh38"         # Reference genome index name
```
How to modify:
- Set thread numbers according to computational resources
- Update `index_name` to match reference genome

**group parameters**
```yaml
  analysis_groups:                # Group analysis settings
    group_names: ["RNA_A", "RNA_N"]  # List of all experimental groups
    sample_replicates: [2, 2]           # Number of replicates per group
```
How to modify:
- Modify group names and replicate numbers to match experimental design
**differential expression genes parameters**
```yaml
  deseq2:                         # Differential expression analysis
    control_group: "RNA_N"      # Control/reference condition
    treatment_group: "RNA_A"    # Treatment/test condition
```
How to modify:
- according experimental design to modify group

#### **Resources Section**
Specify computational resource requirements:

```yaml
resources:
  hisat2_threads: 5      # CPU threads for alignment
  memory_gb: 100         # Maximum memory usage (GB)
```
How to modify: Adjust based on system capabilities and dataset size.

### **Step 3: Validate Configuration**

Before running the pipeline, ensure:

1. All file paths are correct and accessible
2. Sample IDs match FASTQ file prefixes
3. Reference files are properly formatted and indexed
4. Group names are consistent throughout the configuration
5. Resource allocations are within system limits

### **Step 4: Run the Pipeline**

Once configuration is ready, execute:

```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory
snakemake --snakefile ./scripts/snakemake/Preprocess_RNAseq_data/RNAseq_snakefile.py --configfile ./scripts/snakemake/Preprocess_RNAseq_data/config.yaml -j 10

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
snakemake --snakefile ./scripts/snakemake/Preprocess_RNAseq_data/RNAseq_snakefile.py --configfile ./scripts/snakemake/Preprocess_RNAseq_data/config.yaml -j 10
```

## **List of processes**
- [RNAseq_quality_control](#rnaseq_quality_control)
  - [Trim_galore](#trim_galore)
  - [Multiqc](#multiqc)
- [RANseq_mapping_reads_to_genome](#ranseq-mapping-reads-to-genome)
  - [Align_RNA-seq_reads_to_the_reference](#align-rna-seq-reads-to-the-reference-genome-using-hisat2)
  - [Filter_multi-mapped_reads](#filter-multi-mapped-reads)
  - [Convert_sam_files_to_bam_format](#convert-sam-files-to-bam-format-using-samtools)
  - [Sort_bam_files_by_genome_coordinates](#sort-bam-files-by-genome-coordinates)
  - [Mark_and_remove_duplicate_reads](#mark-and-remove-duplicate-reads-using-picard)
  - [Index_bam_files](#index-bam-files-using-samtools)
- [RSeQC_gene_body_coverage](#rseqc_gene_body_coverage)
- [RNAseq_alignment_rates](#rnaseq_alignment_rates)
- [StringTie_expression_quantification](#stringtie_expression_quantification)
- [RANseq_extract_expression_data](#ranseq_extract_expression_data)
- [RNAseq_cluster_analysis](#rnaseq_cluster_analysis)
- [RANseq_PCA_analysis](#ranseq_pca_analysis)
- [RNAseq_differential_expression_analysis](#rnaseq_differential_expression_analysis)
  - [Construct_samplelist](#construct_samplelist)
  - [Prepare_differential_expression_analysis](#prepare_differential_expression_analysis)
  - [Differential_expression_analysis](#differential_expression_analysis)

## **Preprocessing**

## RNAseq_quality_control

### **Trim_galore**

**Description**

Trim_galore is used to remove low quality bases and adapter sequences from the raw RNA-seq reads. It performs quality control by trimming low-quality ends from reads in addition to adapter removal.

**Inputs**
- Raw FASTQ files: `*_1_fq`, `*_2.fq` (for paired-end)
- Raw FASTQ files: `*_1_.fq` (for single-end)

**Parameters**
- `--quality`: Quality threshold for trimming low-quality ends (default: 25)
- `--paired33`:  Indicate that the input is paired-end data with Phred+33 quality score encoding.
- `--stringency`: Overlap with adapter sequence required to trim (default: 3)
- `--length`: Discard reads shorter than this length (default: 35)
- `--fastqc`: Quality Control Report

**Outputs**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)
- FastqQC reports: `*fastqc.html`

**Output directory** 
- `./datasets/RNAseq/quality_control_file/`

### **Multiqc**
**Description**

MultiQC aggregates results from various bioinformatics analyses into a single report. It collects the output from Trim_galore and other QC tools to generate comprehensive quality control metrics.

**Input**
- FastqQC reports: `*fastqc.html`

**Outputs**  
- HTML report: `multiqc_report.html`
- Data directory: `multiqc_data/`

**Output directory**  
- `./datasets/RNAseq/quality_control_file/multiqc_report.html`

## RANseq mapping reads to genome

### **Align RNA-seq reads to the reference genome using HISAT2**

**Description**  
HISAT2 is used to align RNA-seq reads to a reference genome. It is designed specifically for RNA-seq data and can handle spliced alignments efficiently.

**Input**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)

**Parameters**  
- `-x`: Path to HISAT2 index
- `-1`: Path to forward reads
- `-2`: Path to reverse reads (for paired-end)
- `-S`: Output SAM file
- `--dta`: Output alignments tailored for transcript assemblers
- `-p`: Number of threads
- `--summary-file`: Output alignments rate

**Outputs**  
- SAM files: `*accepted_hits.sam`
- Alignment summary: `mapping_summary.txt`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Filter multi-mapped reads**

**Description**  
Use "grep" to process the alignment results and filter out reads that align to multiple genomic locations,retaining only uniquely mapped reads to improve downstream analysis accuracy.

**Parameters**  
- `-v`: nvert match — select lines that do not match the pattern.
- `-E`: Use extended regular expressions (ERE) for the pattern.
- `-w`: Match only whole words (pattern must match an entire word).

**Outputs**  
- Filtered BAM files: `*accepted_hits_NHi1.sam`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Convert sam files to bam format using samtools**

**Description**

Samtools view converts SAM files to BAM format, which is a compressed binary version of the SAM format. This step reduces file size and prepares the data for downstream processing. It can also be used to filter reads based on various criteria.

**Input**
- SAM files: `*accepted_hits_NHi1.sam` 

**Parameters**
- `-S`: Input is SAM format
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**
- BAM files: `*accepted_hits_NHi1.bam`

**Output directory**
- `./datasets/RNAseq/hisat2file/`

### **Sort bam files by genome coordinates**

**Description**  
This step sorts the BAM files by genomic coordinates, which is required for many downstream analyses and improves the efficiency of data access.

**Input**
- BAM files: `*accepted_hits_NHi1.bam`

**Parameters**  
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**  
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam `

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Mark and remove duplicate reads using Picard**

**Description**  
Picard's MarkDuplicates identifies and flags or removes duplicate reads that may have resulted from PCR amplification during library preparation.

**INput**
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam `

**Parameters**  
- `-Xmx15g`: Set the maximum Java heap memory to 15 GB.
- `I`: Input BAM file
- `O`: Output BAM file
- `METRICS_FILE`: File to write duplicate metrics
- `REMOVE_DUPLICATES`: Whether to remove (true) or just mark (false) duplicates(default: true)
- `ASSUME_SORT_ORDER`: Specify the expected sort order of the input BAM file, skipping sort-order validation (default:coordinate)

**Outputs**  
- Deduplicated BAM files: `*[sampleID].bam`
- Duplicate metrics: `*.metricsFile`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

### **Index bam files using samtools**

**Description**  
This step creates an index for the BAM files, which allows for quick random access to the data, essential for many downstream analysis tools.
**Input**
- Deduplicated BAM files: `*[sampleID].bam`

**Outputs**  
- BAM index files: `*[sampleID].bam.bai`

**Output directory**  
- `./datasets/RNAseq/hisat2file/`

## RSeQC_gene_body_coverage

**Description**  
RSeQC's gene body coverage analysis evaluates the uniformity of coverage across gene bodies, which can help identify potential biases in the RNA-seq experiment.

**Input**
- All samples BAM files: `*[sampleID].bam`

**Parameters**  
- `-i`: Input BAM file
- `-r`: Reference gene model in BED format
- `-o`: Output prefix

**Outputs**  
- Gene body coverage plot: `*_geneBodyCoverage.curves.pdf`
- Gene body heatMap plot: `*_geneBodyCoverage.heatMap.pdf`

**Output directory**  
- `./datasets/RNAseq/`

## RNAseq_alignment_rates

**Description**  
This step will present the alignment rates of all samples in a stacked bar chart, providing a quick overview of the mapping quality.

**Input**
- Alignment summary: `mapping_summary.txt` for all samples

**Outputs**  
- Alignment rate plot: `RNAseq_alignment_rates_plot.pdf`

**Output directory**  
- `/datasets/RNAseq/hisat2file`

## StringTie_expression_quantification

**Description**
StringTie is used to estimate transcript abundances from RNA-seq alignments. It performs reference-guided quantification to generate gene and transcript expression levels for downstream analysis.

**Input**
All samples BAM files: `*[sampleID].bam`

**Parameters**

- `-p`: Number of processing threads (e.g., 5)
- `-e`: Quantification mode only (no novel transcript assembly)
- `-B`: Enable output for Ballgown downstream analysis
- `-G`: Reference annotation file in GTF format
- `-A`: Output file for gene abundances in tab-delimited format
- `-o`: Output file for assembled transcripts in GTF format

**Outputs**

- Gene abundance table: `gene_abund.tab`
- Assembled transcripts: `transcripts.gtf`

**Output directory**
`./datasets/RNAseq/stringtiefile/`

## RANseq_extract_expression_data

**Description**  
In this step, we use a custom script to extract the gene expression data (TPM) from the results of Stringtie.

**Input**
- Gene abundance table: `gene_abund.tab`

**Outputs**  
- Expression matrix: `gene_TPM.txt`


**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RNAseq_cluster_analysis

**Description**  
This step performs hierarchical clustering of samples based on their gene expression matrix to identify similarities and differences between samples.

**Input**
- Expression matrix: `gene_TPM.txt`

**Outputs**  
- Clustering dendrogram: `RNAseq_sample_clustering.pdf`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RANseq_PCA_analysis

**Description**  
Principal Component Analysis (PCA) is performed to visualize the overall structure of the data.We employ a custom code that will produce both 2D and 3D PCA maps.

**Input**  
- Expression matrix: `gene_TPM.txt`

**Outputs**  
- 2D PCA plot: `2DPCA_PC1_PC2.pdf` `2DPCA_PC1_PC3.pdf` `2DPCA_PC2_PC3.pdf`
- 2D PCA plot: `3DPCA.pdf.pdf`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

## RNAseq_differential_expression_analysis

### **Construct_samplelist**

**Description**  
DEseq2 requires defining sample names and experimental conditions in advance. This step creates a sample list file that defines the experimental design and the input file path.

**Input**
- Assembled transcripts: `transcripts.gtf`

**Outputs**  
- Sample list file: `samplelist.txt`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

### **Prepare_differential_expression_analysis**

**Description**  
This step prepares the data for differential expression analysis by normalizing the counts and creating the necessary data structures.

**Inputs**  
- Sample list file: `samplelist.txt`

**Outputs**  
- Normalized counts: `gene_count_matrix.count`

**Output directory**  
- `/datasets/RNAseq/stringtiefile`

### **Differential_expression_analysis**

**Description**  
This step identifies differentially expressed genes between conditions using statistical methods DESeq2 

**Inputs**  
- Normalized counts: `gene_count_matrix.count`

**Outputs**  
- Differential expression results: `Exper_vs_ctr_DEG.csv.txt`

**Output directory**  

- `/datasets/RNAseq/stringtiefile`
