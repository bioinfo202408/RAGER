# **RAGER table of contents**
1. [Quick start](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md)
3. [Preprocess ATACseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_ATACseq/ATACseq_analysis.md) 
4. [Joint analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Joint_analysis/Joint_analysis.md)
5. [Custom analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Custom_analysis/Custom_analysis.md)

## **Quick Start**

Welcome to RAGER! This section will guide through the essential steps to configure and run ATAC-seq analysis pipeline.

### **Step 1: Open the Configuration File**

Open the `config.yaml` file in the species-specific directory:

```bash
cd ~/PROJECT/RAGER  #Change the directory to PROJECT/RAGER as the current directory
#mouse
vim ./mouse/scripts/snakemake/Preprocess_ATACseq_data/config.yaml  
#human
vim ./human/scripts/snakemake/Preprocess_ATACseq_data/config.yaml 
# you can use any text editor 
```

### **Step 2: Modifying the Configuration File**

The `config.yaml` file contains all parameters needed to customize RAGER for specific ATAC-seq analysis. Below is a detailed explanation of each section:

#### **Paths Section**
This section defines the directory structure for ATAC-seq analysis:

```yaml
paths:
  base_dir: "datasets"                        # Main data directory
  scripts_dir: "scripts"                      # Directory containing analysis scripts
  atacseq_base_dir: "datasets/ATACseq"        # ATAC-seq analysis base directory
  fastq_dir: "datasets/ATACseq/fastqfile"     # Raw FASTQ files location
  qc_dir: "datasets/ATACseq/quality_control_file"  # Quality control outputs
  bowtie2_dir: "datasets/ATACseq/bowtie2file"      # Alignment outputs
  macs2_dir: "datasets/ATACseq/macs2file"          # Peak calling outputs
  index_dir: "reference/bowtie2index"              # Reference genome index
  geneanno_dir: "reference/geneanno"               # Gene annotation directory
```

**How to modify**: All directories exist or will be created by the pipeline, it is recommended to use the default path above.

#### **References Section**
Specify the reference files for target organism:

```yaml
#mouse
references:
  gtf_file: "reference/geneanno/gencode.vM25.annotation.gtf"  # Gene annotation file
  bed_file: "reference/geneanno/gencode.vM25.annotation.bed"  # Gene annotation in BED format
  bowtie2_index: "GRCm38"                                     # Bowtie2 index name
#human
references:
  gtf_file: "reference/geneanno/gencode.v44.annotation.gtf"   # Gene annotation file
  bed_file: "reference/geneanno/gencode.v44.annotation.bed"   # Gene annotation in BED format
  bowtie2_index: "GRCh38"                                     # Bowtie2 index name
```

**How to modify**: 
- Replace with species-specific GTF and BED files
- Update bowtie2_index to match your reference genome

#### **Samples Section**
Define experimental samples and their properties:

```yaml
samples:
  SRR4032269:           # Sample ID (must match FASTQ file prefix)
    type: "paired"       # "paired" or "single" end sequencing
    group: "ATAC_A"     # Experimental group/condition
    label: "ATAC_A1"    # Human-readable sample label
  SRR4032270: 
    type: "paired"
    group: "ATAC_A"
    label: "ATAC_A2"
  SRR4032271: 
    type: "paired"
    group: "ATAC_N"
    label: "ATAC_N1"
  SRR4032272: 
    type: "paired"
    group: "ATAC_N"
    label: "ATAC_N2"
```

**How to modify**:
- Replace sample IDs with actual sample names
- Set `type` to `"paired"` for paired-end or `"single"` for single-end data
- Assign meaningful `group` names for experimental conditions
- Provide descriptive `label` names for easier interpretation

#### **Parameters Section**

**Trim Galore Parameters:**
```yaml
params:
  trim_galore:                    # Quality control parameters
    phred: 33                     # Quality score encoding
    quality: 25                   # Quality threshold for trimming
    length_paired: 35             # Minimum read length (paired-end)
    length_single: 30             # Minimum read length (single-end)
    stringency_paired: 3          # Adapter stringency (paired-end)
    stringency_single: 1          # Adapter stringency (single-end)
```

**How to modify**: Adjust quality parameters based on data quality requirements. For details, please refer to [trim galore website](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

**Alignment Parameters:**
```yaml
#mouse
  bowtie2:                        # Alignment parameters
    threads: 5                    # Number of CPU threads
    index_name: "GRCm38"         # Reference genome index name
#human
  bowtie2:                        # Alignment parameters
    threads: 5                    # Number of CPU threads
    index_name: "GRCh38"         # Reference genome index name
```

**How to modify**: Set thread numbers according to computational resources and update `index_name` to match reference genome

**Coverage Analysis Parameters:**
```yaml
  bamcoverage:                    # Coverage analysis parameters
    bin_size: 10                  # Bin size for coverage calculation
    processors: 5                 # Number of processors
    normalization: "BPM"          # Normalization method (BPM, RPKM, etc.)
```

**How to modify**: 
- Adjust `bin_size` for coverage resolution (smaller values = higher resolution)
- Set `processors` according to computational resources
- Choose normalization method: "BPM" (bins per million), "RPKM", "RPGC", "None"

It is recommended to use the default parameters above.For details, please refer to [deeptools](https://deeptools.readthedocs.io/en/latest/)

**Calculate access signal Parameters:**
```yaml
  computematrix:                  # Matrix computation parameters
    before: 1000                  # Bases before reference point
    after: 1000                   # Bases after reference point
    bin_size: 10                  # Bin size for matrix
    processors: 10                # Number of processors
```

**How to modify**: 
- Adjust `before` and `after` values to define the region around reference points
- Set `bin_size` for matrix resolution
- Configure `processors` based on available computational resources

It is recommended to use the default parameters above.For details, please refer to [deeptools](https://deeptools.readthedocs.io/en/latest/)

**Peak Calling Parameters:**
```yaml
  macs2:                          # Peak calling parameters
    genome_size: "mm"            # Genome size ("mm" for mouse, "hs" for human)
    qvalue: 0.05                 # Q-value threshold
    keep_dup: "all"              # Duplicate read handling
```

**How to modify**: 
- Change `genome_size` to `"hs"` for human samples, `"mm"` for mouse
- Adjust `qvalue` threshold for peak significance (lower = more stringent)
- Set `keep_dup` to control duplicate reads: "all", "auto", or integer value
For details, please refer to [MACS2](https://pypi.org/project/MACS2/)

**Group Analysis Parameters:**
```yaml
  analysis_groups:  # Group analysis settings
    group_names: ["ATAC_A", "ATAC_N"] # List of all experimental groups
    sample_labels: ["ATAC_A1", "ATAC_A2", "ATAC_N1", "ATAC_N2"] # Sample labels
    column_ranges: ["1:200", "201:400", "401:600", "601:800"] # Column ranges
    sample_replicates: [2, 2] # Number of replicates per group
    sample_order: [1, 2, 3, 4]  # Sample ordering for visualization
```

**How to modify**: 
- Modify `group_names` and `sample_labels` to match experimental design
- Adjust `column_ranges` for the number of samples. One more sample increases by 200
- Update `sample_replicates` to reflect actual replicate numbers
- Set `sample_order` to control sample arrangement in visualizations

**Differential Peaks Analysis Parameters:**
```yaml
  differential_analysis:          # Differential accessibility analysis
    treatment_group: "ATAC_A"    # Treatment condition
    control_group: "ATAC_N"      # Control/reference condition
    treatment_samples: ["SRR4032269", "SRR4032270"]  # Treatment sample IDs
    control_samples: ["SRR4032271", "SRR4032272"]    # Control sample IDs
```

**How to modify**: Adjust according to experimental design and sample IDs

#### **Output Directories Section**
Define output subdirectory structure:

```yaml
output_dirs:
  bowtie2_subdir: "bowtie2file"   # Subdirectory for alignment outputs
  macs2_subdir: "macs2file"       # Subdirectory for peak calling outputs
```

**How to modify**: Change subdirectory names if needed, it is recommended to use the default parameters.

### **Step 3: Validate Configuration**

Before running the pipeline, ensure:

1. All file paths are correct and accessible
2. Sample IDs match FASTQ file prefixes
3. Bowtie2 reference index is properly built and indexed
4. Group names are consistent throughout the configuration
5. Resource allocations are within system limits
6. Treatment and control samples are correctly specified
7. Column ranges and sample order match your experimental design

### **Step 4: Run the Pipeline**

Once configuration is ready, execute:

```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory
snakemake --snakefile ./scripts/snakemake/Preprocess_ATACseq_data/ATACseq_snakefile.py --configfile ./scripts/snakemake/Preprocess_ATACseq_data/config.yaml -j 10

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
snakemake --snakefile ./scripts/snakemake/Preprocess_ATACseq_data/ATACseq_snakefile.py --configfile ./scripts/snakemake/Preprocess_ATACseq_data/config.yaml -j 10
```
## **List of processes**
- [ATAC_quality_control](#)
  - [Trim_galore](#)
  - [Multiqc](#)
- [ATACseq_mapping_reads_to_genome](#)
  - [Align_ATACseq_reads_to_the_reference](#)
  - [Filter_multi-mapped_reads](#)
  - [Convert_sam_to_bam_format](#)
  - [Sort_bam_files_by_genome_coordinates](#)
  - [Mark_and_remove_duplicate_reads](#)
  - [Index_bam_files](#)
- [ATACseq_alignment_rates](#)
- [ATACseqQC](#)
- [Trans_gene_anno_to_bed](#)
- [Bam_to_bw](#)
- [Quantify TPM signals](#)
- [ATACseq_cluster_analysis](#)
- [ATACseq_PCA_analysis](#)
- [Merge_bam_file](#)
- [ATACseq_peak_calling](#)
- [ATACseq_peak_annotation](#)
- [ATACseq_peak_bar_plot](#)

## **Preprocessing**

## ATAC_quality_control

### **Trim_galore**

**Description**

Trim_galore is used to remove low quality bases and adapter sequences from the raw ATAC-seq reads. It performs quality control by trimming low-quality ends from reads in addition to adapter removal, which is particularly important for ATAC-seq data to ensure accurate mapping and peak calling.

**Inputs**
- Raw FASTQ files: `*_1.fq`, `*_2.fq` (for paired-end)
- Raw FASTQ files: `*_1.fq` (for single-end)

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
- `./datasets/ATACseq/quality_control_file/`

### **Multiqc**

**Description**

MultiQC aggregates results from various bioinformatics analyses into a single report. It collects the output from Trim_galore and other QC tools to generate comprehensive quality control metrics for ATAC-seq data.

**Input**
- FastqQC reports: `*fastqc.html`

**Outputs**  
- HTML report: `multiqc_report.html`
- Data directory: `multiqc_data/`

**Output directory**  
- `./datasets/ATACseq/quality_control_file/`

## ATACseq_mapping_reads_to_genome

### **Align_ATACseq_reads_to_the_reference**

**Description**  
Bowtie2 is used to align ATAC-seq reads to a reference genome. It is well-suited for short read alignment and is commonly used for ATAC-seq data to map open chromatin regions.

**Input**
- Trimmed FASTQ files: `*_1_val_1.fq`, `*_2_val_2.fq` (for paired-end)
- Trimmed FASTQ files: `*_1_trimmed.fq` (for single-end)

**Parameters**  
- `-x`: Reference genome index prefix
- `-1`: Forward/left reads file (R1)
- `-2`: Reverse/right reads file (R2)
- `-t`: Print wall-clock time taken by search phases
- `-q`: Input reads are in FASTQ format
- `-N`: Max mismatches in seed alignment (default: 1)
- `-L`: Length of seed substrings (default: 25)
- `-S`: Output alignment results in SAM format
- `-p`: Launch specified number of parallel search threads
- `--no-mixed`: Suppress unpaired alignments for paired reads
- `--no-discordant`: Suppress discordant alignments for paired reads

**Outputs**  
- SAM files: `*accepted_hits.sam`
- Alignment summary: `mapping_summary.txt`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Filter_multi-mapped_reads**

**Description**  
This step filters out reads that align to multiple genomic locations, retaining only uniquely mapped reads to improve the accuracy of peak calling and other downstream analyses.

**Parameters**  
- `-H`: Extract header lines only
- `grep 'AS:'`: Find reads with alignment score (AS:) tag
- `grep -v 'XS:'`: Exclude reads with suboptimal alignment score (XS:) tag
**Outputs**  
- Filtered SAM files: `*accepted_hits_NHi1.sam`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Convert_sam_to_bam_format**

**Description**

Samtools view converts SAM files to BAM format, which is a compressed binary version of the SAM format. This step reduces file size and prepares the data for downstream processing.

**Input**
- Filtered SAM files: `*accepted_hits_NHi1.sam`

**Parameters**
- `-S`: Input is SAM format
- `-b`: Output BAM format
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**
- BAM files: `*accepted_hits_NHi1.bam`

**Output directory**
- `./datasets/ATACseq/bowtie2file/`

### **Sort.bam_files_by_genome_coordinates**

**Description**  
This step sorts the BAM files by genomic coordinates, which is required for many downstream analyses including peak calling and visualization.

**Input**
- BAM files: `*accepted_hits_NHi1.bam`

**Parameters**  
- `-o`: Output file name
- `-@`: Number of threads

**Outputs**  
- Sorted BAM files: `*accepted_hits_NHi1.sorted.bam`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Mark_and_remove_duplicate_reads_using_Picard**

**Description**  
Picard's MarkDuplicates identifies and flags or removes duplicate reads that may have resulted from PCR amplification during library preparation, which is particularly important for ATAC-seq data to prevent bias in peak calling.

**Input**
- Sorted BAM files: `*accepted_hits_NHi1_sorted.bam`

**Parameters**  
- `-Xmx15g`: Set the maximum Java heap memory to 15 GB
- `I`: Input BAM file
- `O`: Output BAM file
- `METRICS_FILE`: File to write duplicate metrics
- `VALIDATION_STRINGENCY`: Controls validation strictness level(default: LENIENT)
- `REMOVE_DUPLICATES`: Whether to remove (true) or just mark (false) duplicates(default: true)
- `ASSUME_SORT_ORDER`: Specify the expected sort order of the input BAM file, skipping sort-order validation (default:coordinate)

**Outputs**  
- Deduplicated BAM files: `*[sampleID].bam`
- Duplicate metrics: `*.metricsFile`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

### **Index_bam_files_using_samtools**

**Description**  
This step creates an index for the BAM files, which allows for quick random access to the data, essential for visualization and downstream analyses.

**Input**
- Deduplicated BAM files: `*[sampleID].bam`

**Outputs**  
- BAM index files: `*[sampleID].bam.bai`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_alignment_rates

**Description**  
This step will present the alignment rates of all samples in a stacked bar chart, providing a quick overview of the mapping quality for ATAC-seq data.

**Input**
- Alignment summary: `mapping_summary.txt` for all samples

**Outputs**  
- Alignment rate plot: `ATACseq_alignment_rates_plot.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseqQC

**Description**  
ATACseqQC performs quality control analyses specific to ATAC-seq data, including fragment size distribution, TSS enrichment, and heatmap.

**Input**
- Indexed BAM files: `*[sampleID].bam` and `*[sampleID].bam.bai`

**Outputs**  
- Fragment size distribution plot: `fragmentSizeDistribution.pdf`
- Heatmap of signal around genomic features: `heatmap.pdf`
- Coverage curve plot: `coverage_curve_plot.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Trans_gene_anno_to_bed

**Description**  
This step converts gene annotation from GTF format to BED format, which is required for many downstream analyses including peak annotation.

**Input**
- Gene annotation file: `gencode.v44.annotation.gtf`

**Outputs**  
- BED format gene annotation: `gencode.v44.annotation.bed`

**Output directory**  
- `./reference_annotation_file/geneanno/`

## Bam_to_bw

**Description**  
This step uses the **Deeptool** software bamCoverage tool to convert BAM files into bigWig format, which is more efficient for subsequent calculation of TPM signals and visualization of the genome browser.

**Input**
- Indexed BAM files: `*[sampleID].bam` and `*[sampleID].bam.bai`

**Parameters**
- `--binsize`: Bin size for signal calculation (default:10)
- `--ignoreDuplicates`: Skip PCR duplicate reads marked by Picard
- `--normalizeUsing` : Coverage normalization strategy (default:BPM)
- `-p`: Number of threads
- `-of`: Output format (bigwig)

**Outputs**  
- BigWig files: `*[sampleID].bw`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Quantify TPM signals

**Description**  
This step uses the **Deeptool** software computeMatrix tool to quantifies the normalized signal (TPM - Tags Per Million) from ATAC-seq data, which allows for comparison between samples.

**Input**
- Indexed BAM files: `*[sampleID].bam`

**Parameters**
- `--reference-point`: Reference-point centered analysis mode
- `-b`: Distance upstream of reference point (bp) (default:1000)
- `-a`: Distance downstream of reference point (bp) (default:1000)
- `--binsize`: Bin size for signal calculation (default:10)

**Outputs**  
- Normalized signal table: `all_scaled.tab`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_cluster_analysis

**Description**  
This step performs hierarchical clustering of samples based on their ATAC-seq signals to identify similarities and differences between samples.

**Input**
- Normalized signal table: `all_scaled.tab`

**Outputs**  
- Clustering dendrogram: `sample_clustering.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_PCA_analysis

**Description**  
Principal Component Analysis (PCA) is performed to visualize the overall structure of the ATAC-seq data and identify major sources of variation.We employ a custom code that will produce both 2D and 3D PCA maps.

**Input**  
- Normalized signal table: `all_scaled.tab`

**Outputs**  
- 2D PCA plot: `2DPCA_PC1_PC2.pdf` `2DPCA_PC1_PC3.pdf` `2DPCA_PC2_PC3.pdf`
- 2D PCA plot: `3DPCA.pdf.pdf`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## Merge_bam_file

**Description**  
This step merges BAM files from replicates or related samples to increase the signal-to-noise ratio for peak calling.

**Input**
- Indexed BAM files: `*[sampleID].bam`

**Outputs**  
- Merged BAM files: `merged.bam`
- Merged BAM index files: `merged.bam.bai`

**Output directory**  
- `./datasets/ATACseq/bowtie2file/`

## ATACseq_peak_calling

**Description**  
MACS2 is used to call peaks from ATAC-seq data, which represent regions of open chromatin. This step identifies differential peaks between conditions.

**Input**
- Merged BAM files: `experiment.bam`, `ctr.bam`

**Parameters**
- `--treatment`: Treatment BAM file
- `--control`: Control BAM file
- `--genome`: Species of this data
- `--name`: Output file prefix
- `--outdir`: Output directory
- `-q`: Q-value cutoff for peak detection

**Outputs**  
- Peak files: `ATACexperiment_vs_ATACctr_summits.bed`, `ATACctr_vs_ATACexperiment_summits.bed`

**Output directory**  
- `./datasets/ATACseq/macs2file/`

## ATACseq_peak_annotation

**Description**  
This step annotates the called peaks with genomic features like promoters, enhancers, and gene bodies to provide functional context.

**Input**
- Peak files: `ATACexperiment_vs_ATACctr_summits.bed`, `ATACctr_vs_ATACexperiment_summits.bed`

**Outputs**  
- Annotated peak tables: `Experiment_vs_ctr_peak.csv`, `Ctr_vs_experiment_peak.csv`

**Output directory**  
- `./datasets/ATACseq/macs2file/`

## ATACseq_peak_bar_plot

**Description**  
This step creates bar plots to visualize the number and distribution of differential peaks between conditions.

**Input**
- Annotated peak tables: `GDC_vs_ctr_peak.csv`, `ctr_vs_GDC_peak.csv`

**Outputs**  
- Diverging bar plot: `peaks_diverging_bar.pdf`
- All peaks bar plot: `all_peaks_bar.pdf`

**Output directory**  
- `./datasets/ATACseq/macs2file/`
