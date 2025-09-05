# **RAGER table of contents**
1. [Quick start](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md)
3. [Preprocess ATACseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_ATACseq/ATACseq_analysis.md) 
4. [Joint analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Joint_analysis/Joint_analysis.md)
5. [Custom analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Custom_analysis/Custom_analysis.md)

## **Quick Start**

Welcome to RAGER! This section will guide through the essential steps to configure and run promoter and enhancer region analysis pipeline.

### **Step 1: Open the Configuration File**

Open the `config.yaml` file in the species-specific directory:

```bash
cd ~/PROJECT/RAGER  #Change the directory to PROJECT/RAGER as the current directory
#mouse
vim ./mouse/scripts/snakemake/Joint_analysis/config.yaml  
#human
vim ./human/scripts/snakemake/Joint_analysis/config.yaml 
# you can use any text editor 
```

### **Step 2: Modifying the Configuration File**

The `config.yaml` file contains all parameters needed to customize RAGER for promoter and enhancer region analysis. Below is a detailed explanation of each section:

#### **Paths Section**
This section defines the directory structure for promoter and enhancer analysis:

```yaml
paths:
  base_dir: "datasets"                                    # Main data directory
  scripts_dir: "scripts"                                  # Directory containing analysis scripts
  
  atacseq_dir: "datasets/ATACseq"                        # ATAC-seq data directory
  atacseq_bowtie2_dir: "datasets/ATACseq/bowtie2file"    # ATAC-seq alignment files
  rnaseq_dir: "datasets/RNAseq"                          # RNA-seq data directory
  
  promoter_analysis_dir: "datasets/Promoter_region_analysis"  # Promoter analysis outputs
  enhancer_analysis_dir: "datasets/Enhancer_region_analysis"  # Enhancer analysis outputs
  
  geneanno_dir: "reference"                              # Reference annotation directory
  geneanno_bed: "reference/geneanno/gencode.vM25.annotation.bed"  # Gene annotation BED file
  enhancer_bed: "reference/enhanno/EnhancerAltasv2.0_eRNA.bed"    # Enhancer annotation BED file
  
  jaspar_raw: "datasets/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"  # Raw JASPAR database
  jaspar_dir: "datasets/JASPAR2024_CORE"                 # Processed JASPAR directory
  tf_mapping_file: "datasets/Mus_musculus_TF.txt"        # Transcription factor mapping file
```

**How to modify**: 
- All directories will be created by the pipeline if they don't exist
- Update `geneanno_bed` for different species (e.g., `reference/geneanno/gencode.v44.annotation.bed` for human)
- Update `tf_mapping_file` for different species (e.g., `datasets/Homo_sapiens_TF.txt` for human)
- Except for `geneanno_bed` and `tf_mapping_file`,it is recommended to use the default paths above

#### **Input Files Section**
Specify the input files from previous ATAC-seq and RNA-seq analyses:

```yaml
input_files:
  up_peak_csv: "datasets/ATACseq/macs2file/ATAC_A_vs_ATAC_N_peak.csv"              # Upregulated peaks annotation
  down_peak_csv: "datasets/ATACseq/macs2file/ATAC_N_vs_ATAC_A_peak.csv"            # Downregulated peaks annotation
  
  up_peak_xls: "datasets/ATACseq/macs2file/ATAC_A_vs_ATAC_N/ATAC_A_vs_ATAC_N_peaks.xls"      # Upregulated peak details
  down_peak_xls: "datasets/ATACseq/macs2file/ATAC_N_vs_ATAC_A/ATAC_N_vs_ATAC_A_peaks.xls"    # Downregulated peak details
  
  deg_csv: "datasets/RNAseq/stringtiefile/A_vs_N_DEG.csv"                          # Differential expressed genes
```

**How to modify**: 
- Update file paths to match your actual ATAC-seq and RNA-seq analysis output locations
- Ensure the peak files contain the required differential accessibility information
- Verify that DEG file contains differential expression analysis results

#### **Samples Section**
Define experimental samples and their properties:

```yaml
samples:
  SRR4032269:           # Sample ID
    group: "ATAC_A"     # Experimental group/condition
    label: "ATACA1"     # Human-readable sample label
  SRR4032270:
    group: "ATAC_A"
    label: "ATACA2"
  SRR4032271:
    group: "ATAC_N"
    label: "ATACN1"
  SRR4032272:
    group: "ATAC_N"
    label: "ATACN2"
```

**How to modify**:
- Replace sample IDs with actual sample names from your ATAC-seq experiment
- Assign meaningful `group` names for experimental conditions
- Provide descriptive `label` names for easier interpretation

#### **Parameters Section**

**Differential Expression Thresholds:**
```yaml
params:
  log2fc_up: 1                    # Log2 fold change threshold for upregulation
  log2fc_down: -1                 # Log2 fold change threshold for downregulation
  padj: 0.05                      # Adjusted p-value threshold
```

**How to modify**: 
- Adjust `log2fc_up` and `log2fc_down` for more or less stringent fold change requirements
- Modify `padj` threshold for statistical significance (lower = more stringent)

**Matrix Computation Parameters:**
```yaml
  computematrix:                  # Matrix computation parameters for promoter/enhancer regions
    before: 1000                  # Bases before reference point
    after: 1000                   # Bases after reference point
    bin_size: 10                  # Bin size for matrix
    processors: 5                 # Number of processors
```

**How to modify**: 
- Adjust `before` and `after` values to define the region around promoter/enhancer elements
- Set `bin_size` for matrix resolution (smaller values = higher resolution)
- Configure `processors` based on available computational resources

It is recommended to use the default parameters above.For details, please refer to [deeptools](https://deeptools.readthedocs.io/en/latest/)

**Group Analysis Parameters:**
```yaml
  group_analysis:                 # Group analysis settings
    column_ranges: ["1:200", "201:400", "401:600", "601:800"]    # Column ranges for heatmap
    sample_labels: ["ATACA1", "ATACA2", "ATACN1", "ATACN2"]      # Sample labels
    treatment_group: ["ATACA1", "ATACA2"]                        # Treatment sample labels
    control_group: ["ATACN1", "ATACN2"]                          # Control sample labels
```

**How to modify**: 
- Update `sample_labels` to match the labels defined in the samples section
- Adjust `column_ranges` based on the number of samples (add 200 for each additional sample)
- Modify `treatment_group` and `control_group` to reflect experimental design

**Enrichment Analysis Directories:**
```yaml
  enrichment_analysis:            # Enrichment analysis output directories
    promoter_up_dir: "datasets/Promoter_region_analysis/Enrichment_analysis/Up"      # Promoter upregulated enrichment
    promoter_down_dir: "datasets/Promoter_region_analysis/Enrichment_analysis/Down"  # Promoter downregulated enrichment
    enhancer_up_dir: "datasets/Enhancer_region_analysis/Enrichment_analysis/Up"      # Enhancer upregulated enrichment
    enhancer_down_dir: "datasets/Enhancer_region_analysis/Enrichment_analysis/Down"  # Enhancer downregulated enrichment
```

**How to modify**: Change directory names if needed, or keep defaults for consistency

**Regulation Analysis Directories:**
```yaml
  regulation_analysis:            # Regulation analysis output directories
    promoter_up_dir: "datasets/Promoter_region_analysis/Regulation_analysis/Up"      # Promoter upregulated regulation
    promoter_down_dir: "datasets/Promoter_region_analysis/Regulation_analysis/Down"  # Promoter downregulated regulation
    enhancer_up_dir: "datasets/Enhancer_region_analysis/Regulation_analysis/Up"      # Enhancer upregulated regulation
    enhancer_down_dir: "datasets/Enhancer_region_analysis/Regulation_analysis/Down"  # Enhancer downregulated regulation
```

**How to modify**: Change directory names if needed, or keep defaults for consistency

### **Step 3: Validate Configuration**

Before running the pipeline, ensure:

1. All input files from ATAC-seq and RNA-seq analyses are present and accessible
2. Sample IDs and labels are consistent with previous analyses
3. Reference annotation files (gene and enhancer BED files) are properly formatted
4. Differential expression thresholds are appropriate for your analysis
5. Resource allocations are within system limits
6. Treatment and control groups are correctly specified

### **Step 4: Run the Pipeline**

Once configuration is ready, execute:

```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory
snakemake --snakefile ./scripts/snakemake/Joint_analysis/Analysis_snakefile.py --configfile ./scripts/snakemake/Joint_analysis/config.yaml -j 10

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
snakemake --snakefile ./scripts/snakemake/Joint_analysis/Analysis_snakefile.py --configfile ./scripts/snakemake/Joint_analysis/config.yaml -j 10
```

## **List of processes**
- [Analysis_of_the_promoter_region](#)
    - [Get_shared_up_gene_InPromoter](#)
    - [Get_shared_down_gene_InPromoter](#)
    - [Expression_shared_peak_annotation](#)
    - [Calculate_shared_peak_TPM_signal](#)
    - [Plot_RNAseq_ATACseq_correlation](#)
    - [Enrichment_analysis_in_promoter_region](#)
        - [Enrichment_analysis_for_shared_gene](#)
    - [Regulatory_Analysis_In_promoter_region](#)
        - [Format_TF_motifs](#)
        - [Extract_shared_Gene_promoterSeq](#)
        - [Motif_enrich](#)
        - [Constructe_TF_Gene_network](#)
        - [Extract_enriched_motifs](#)
        - [Plot_TF_heatmap]()
- [Analysis_of_the_enhancer_region]()
    - [Get_RNAseq_DE_gene_annotation](#)
    - [Calculate_enhancer_DEGene_distance](#)
    - [Get_enhancer_annotation](#)
    - [Calculate_enhancer_TPM_signal](#)
    - [Get_shared_Up_Enhancer_and_shared_Up_gene](#)
    - [Get_shared_Down_Enhancer_and_shared_Down_gene](#)
    - [Enrichment_analysis_in_enhancer_region](#)
        - [Enrichment_analysis_for_shared_gene](#)
    - [Regulatory_Analysis_in_enhancer_region](#)
        - [Extract_enhancer_sequences](#)
        - [Motif_enrich](#)
        - [Constructe_TF_Gene_network](#)
        - [Extract_enriched_motifs](#)
        - [Plot_TF_heatmap](#)
## **Preprocessing**
## **Analysis of the promoter region**

### **Get_shared_up_gene_InPromoter**

**Description**

This process identifies genes that are significantly upregulated in RNA-seq data and have accessible chromatin regions (peaks)  from ATAC-seq data. It performs intersection analysis to find overlapping upregulated genes between RNA-seq differential expression results and ATAC-seq peak annotations.

**Inputs**
- `stringtiefile/experiment_vs_ctr_DEG.csv`:RNA-seq differential expression results
- `macs2file/experiment_vs_ctr_peak.csv`:ATAC-seq peak file


**Outputs**
- `UpVenn.pdf`: Venn diagram showing overlap between upregulated genes and accessible promoters
- `shared_up_genes_INpromote.txt`: List of shared upregulated genes in promoter regions

**Output directory**
- `./datasets/Promoter_region_analysis/`

### **Get_shared_down_gene_InPromoter**

**Description**

Similar to the upregulated gene analysis, this process identifies genes that are significantly downregulated in RNA-seq data and have accessible chromatin regions from ATAC-seq data.

**Inputs**
- `stringtiefile/experiment_vs_ctr_DEG.csv`:RNA-seq differential expression results
- `macs2file/ctr_vs_experiment_peak.csv`:ATAC-seq peak file


**Outputs**
- `DownVenn.pdf`: Venn diagram showing overlap between downregulated genes and accessible promoters
- `shared_down_genes_INpromote.txt`: List of shared downregulated genes in promoter regions

**Output directory**
- `./datasets/Promoter_region_analysis/`

### **Extract_shared_peak_annotation**

**Description**

This step will generate the annotation file of the shared peak to prepare for the subsequent calculation of the TPM signal of the peak

**Inputs**
- `shared_up_genes_INpromote.txt`:List of shared upregulated genes in promoter regions
- `macs2file/ATACexperiment_vs_ATACctr/ATACexperiment_vs_ATACctr_peaks.xls`:ATAC-seq peak file
- `shared_down_genes_INpromote.txt`:List of shared downregulated genes in promoter regions
- `/macs2file/ATACctr_vs_ATACexperiment/ATACctr_vs_ATACexperiment_peaks.xls`:ATAC-seq peak file


**Outputs**
- `shared_up_peak_INpromote.bed`: BED file of accessible peaks in upregulated gene
- `shared_down_peak_INpromote.bed`: BED file of accessible peaks in downregulated gene
- `shared_peak_INpromote.bed`: Combined BED file of all shared accessible peaks

**Output directory**
- `./datasets/Promoter_region_analysis/`

### **Calculate_shared_peak_TPM_signal**

**Description**
In this step, the deeptool software computeMatrix tool was used to calculate the TPM signal of the shared peak

**Inputs**
- `shared_peak_INpromote.bed`
- `*[sampleID].bw`:BigWig files

**Parameters**
- `--reference-point`: Reference-point centered analysis mode
- `-b`: Distance upstream of reference point (bp) (default:1000)
- `-a`: Distance downstream of reference point (bp) (default:1000)
- `--binsize`: Bin size for signal calculation (default:10)

**Outputs**
- `shared_peak_scaled.tab`: Tab-delimited file containing normalized TPM signals for all shared peaks

**Output directory**
- `./datasets/Promoter_region_analysis/`

### **Plot_RNAseq_ATACseq_correlation**

**Description**

Generates correlation plots between RNA-seq gene expression levels and ATAC-seq chromatin accessibility signals in promoter regions to assess the relationship between gene expression and chromatin accessibility.

**Inputs**
- `experiment_vs_ctr_DEG.csv`:RNA-seq differential expression results
- `shared_peak_scaled.tab`:TPM signals for all shared peaks
- `shared_up_genes_INpromote.txt`:List of shared upregulated genes in promoter regions
- `shared_down_genes_INpromote.txt`:List of shared downregulated genes in promoter regions

**Outputs**
- `Expr_ccess_corr_plot.pdf`: Correlation plot between expression and accessibility

**Output directory**
- `./datasets/Promoter_region_analysis/`

## **Enrichment_analysis_for_shared_gene**

**Description**

Performs Gene Set Enrichment Analysis (GSEA) for both upregulated and downregulated genes that have accessible promoter regions, identifying biological processes and pathways significantly enriched in the gene sets.

**Inputs**
- `shared_up_genes_INpromote.txt`
- `shared_down_genes_INpromote.txt`


**Parameters**
- Statistical significance threshold (p-value < 0.05 |NES| > 1)


**Outputs**
- `sharedUpGene_INPromote_gsea_GOenrich.csv`: GO enrichment results for upregulated genes
- `sharedUpGene_INPromote_gsea_GOenrich.rds`: GO enrichment results for upregulated genes(rds format)
- `sharedUpGene_INPromote_gsea_KEGGenrich.csv`: KEGG enrichment results for upregulated genes
- `sharedUpGene_INPromote_gsea_KEGGenrich.csv`: KEGG enrichment results for upregulated genes(rds format)
- `sharedDownGene_INPromote_gsea_GOenrich.csv`: GO enrichment results for downregulated genes
- `sharedDownGene_INPromote_gsea_GOenrich.csv`: GO enrichment results for downregulated genes(rds format)
- `sharedDownGene_INPromote_gsea_KEGGenrich.csv`: KEGG enrichment results for downregulated genes
- `sharedDownGene_INPromote_gsea_KEGGenrich.csv`: KEGG enrichment results for downregulated genes(rds format)

**Note:** The rds format file is for user custom drawing GSEAplot.

**Output directory**
- `./datasets/Promoter_region_analysis/Enrichment_analysis/Up`
- `./datasets/Promoter_region_analysis/Enrichment_analysis/Down`
## **Regulatory_Analysis_In_promoter_region**
### **Format_TF_motifs**

**Description**

Formats transcription factor motif databases (JASPAR2024) into the appropriate format for motif enrichment analysis, converting motif position weight matrices into MEME format.

**Inputs**
- `JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt`:JASPAR2024 motif database

**Outputs**
- `./datasets/JASPAR2024_CORE/*.meme`: Individual MEME format files for each transcription factor motif

**Output directory**
- `./datasets/JASPAR2024_CORE/`

### **Extract_shared_Gene_promoterSeq**

**Description**

Extracts DNA sequences from promoter regions of shared upregulated and downregulated genes for subsequent motif enrichment analysis.

**Inputs**
- `shared_up_genes_INpromote.txt`
- `shared_down_genes_INpromote.txt`
- Reference genome FASTA file

**Parameters**
- Promoter region definition ( TSS upstream=2000 , downstream=20)

**Outputs**
- `Regulation_analysis/Up/promoter_seqs.fa`: FASTA sequences for upregulated gene promoters
- `Regulation_analysis/Down/promoter_seqs.fa`: FASTA sequences for downregulated gene promoters

**Output directory**
- `./datasets/Promoter_region_analysis/Regulation_analysis/`

### **Motif_enrich**

**Description**

We performed transcription factor motif enrichment analysis of promoter sequences using the ame tool of MEME software to identify transcription factors that were significantly enriched in accessible promoter regions of shared genes.

**Inputs**
- `Regulation_analysis/Up/promoter_seqs.fa`
- `Regulation_analysis/Down/promoter_seqs.fa`
- `./datasets/JASPAR2024_CORE/*.meme`:JASPAR2024 motif database in MEME format


**Outputs**
- `Regulation_analysis/Up/motif_enrich_result/`
- `Regulation_analysis/Down/motif_enrich_result`
**Output directory**
- `./datasets/Promoter_region_analysis/Regulation_analysis/`

### **Constructe_TF_Gene_network**

**Description**

Constructs transcription factor-gene regulatory networks based on motif enrichment results , linking transcription factors to their potential target genes.

**Inputs**
- `Regulation_analysis/Up/motif_enrich_result/`: Motif enrichment results
- `Regulation_analysis/Down/motif_enrich_result`:Motif enrichment results


**Outputs**
- `Regulation_analysis/Up/SharedUPgene_in_promoter_TF_gene_network.txt`: TF-gene network for upregulated genes
- `Regulation_analysis/Down/SharedDOWNgene_in_promoter_TF_gene_network.txt`: TF-gene network for downregulated genes

**Output directory**
- `./datasets/Promoter_region_analysis/Regulation_analysis/`

### **Extract_enriched_motifs**

**Description**

Extracts and summarizes the list of significantly enriched transcription factor motifs from the motif enrichment analysis results.

**Inputs**
- `Regulation_analysis/Up/motif_enrich_result/`: Motif enrichment analysis results
- `Regulation_analysis/Down/motif_enrich_result/`:Motif enrichment results

**Parameters**
- `FDR < 0.05`:FDR cutoffs for significance


**Outputs**
- `Regulation_analysis/Up/enriched_tfs_list.txt`: List of enriched TFs for upregulated genes
- `Regulation_analysis/Down/enriched_tfs_list.txt`: List of enriched TFs for downregulated genes

**Output directory**
- `./datasets/Promoter_region_analysis/Regulation_analysis/`

### **Plot_TF_heatmap**

**Description**

Generates heatmap visualizations showing the expression patterns of enriched transcription factors across samples, providing insight into TF activity patterns.

**Inputs**
- `Regulation_analysis/Up/enriched_tfs_list.txt`: List of enriched TFs for upregulated genes
- `Regulation_analysis/Down/enriched_tfs_list.txt`: List of enriched TFs for downregulated genes
- `stringtiefile/experiment_vs_ctr_DEG.csv`:RNA-seq differential expression results

**Parameters**
- `Log2FC > 1`: Log2FC cutoffs for TF differential expression
- `FDR < 0.05`: FDR cutoffs for TF differential expression significance

**Outputs**
- `shared_UpGene_InPromote_tfs_heatmap_RNAseqExpr.pdf`: TF expression heatmap for upregulated gene analysis
- `shared_DownGene_InPromote_tfs_heatmap_RNAseqExpr.pdf`: TF expression heatmap for downregulated gene analysis

**Output directory**
- `./datasets/Promoter_region_analysis/Regulation_analysis/`

## **Analysis of the enhancer region**

### **Get_RNAseq_DE_gene_annotation**

**Description**

Annotates differentially expressed genes from RNA-seq analysis with genomic coordinates and generates BED format files for both upregulated and downregulated genes to enable spatial analysis with enhancer regions.

**Inputs**
- `stringtiefile/experiment_vs_ctr_DEG.csv` :RNA-seq differential expression results
- `gencode.v44.annotation.bed` :Gene annotation files (bed fomat)



**Outputs**
- `UpGene_anno.bed`: BED format annotations for upregulated genes
- `DownGene_anno.bed`: BED format annotations for downregulated genes

**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Calculate_enhancer_DEGene_distance**

**Description**

Calculates the genomic distances between enhancer regions and differentially expressed genes, identifying enhancers within a specified distance (1Mb) of upregulated and downregulated genes.

**Inputs**
- `UpGene_anno.bed`: BED format annotations for upregulated genes
- `DownGene_anno.bed`: BED format annotations for downregulated genes
- `EnhancerAtlasv2.0_eRNA.bed` :Annotation of enhancers

**Parameters**
-  `distance < 1Mb` Maximum distance threshold 


**Outputs**
- `enh_UpGene_within_1mb.txt`: Enhancer-gene pairs for upregulated genes within 1Mb
- `enh_DownGene_within_1mb.txt`: Enhancer-gene pairs for downregulated genes within 1Mb

**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Get_enhancer_annotation**

**Description**

Annotates enhancer regions that are associated with differentially expressed genes,create annotation files for these enhancers.

**Inputs**
- `enh_UpGene_within_1mb.txt`: Enhancer-gene pairs for upregulated genes within 1Mb
- `enh_DownGene_within_1mb.txt`: Enhancer-gene pairs for downregulated genes within 1Mb
- `EnhancerAtlasv2.0_eRNA.bed` :Annotation of enhancers


**Outputs**
- `annotation_of_enh_for_potential_regulate_UpGene.bed`: Enahncers annotation of potentially regulated up-regulated genes
- `annotation_of_enh_for_potential_regulate_DownGene.bed`: Enahncers annotation of potentially regulated down-regulated genes
**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Calculate_enhancer_TPM_signal**

**Description**

Calculates normalized TPM signal intensities for enhancer regions associated with differentially expressed genes, enabling identify Co-upregulated enhancers and co-downregulated enhancers.
**Inputs**
- `annotation_of_enh_for_potential_regulate_UpGene.bed`: Enahncers annotation of potentially regulated up-regulated genes
- `annotation_of_enh_for_potential_regulate_DownGene.bed`: Enahncers annotation of potentially regulated down-regulated genes
- `*[sampleID].bw`:BigWig files

**Parameters**
- `--reference-point`: Reference-point centered analysis mode
- `-b`: Distance upstream of reference point (bp) (default:1000)
- `-a`: Distance downstream of reference point (bp) (default:1000)
- `--binsize`: Bin size for signal calculation (default:10)


**Outputs**
- `enh_expr_signal_for_UPgene.tab`: Tab-delimited signal table for enhancers
- `enh_sortedRegions_for_UPgene.bed`: Sorted enhancer regions
- `enh_expr_signal_for_DOWNgene.tab`: Tab-delimited signal table for enhancers
- `enh_sortedRegions_for_DOWNgene.bed`: Sorted enhancer regions
**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Get_shared_Up_Enhancer_and_shared_Up_gene**

**Description**

Identifies enhancers that show consistent accessibility changes and are associated with upregulated genes across samples, creating lists of shared enhancers and their associated genes.

**Inputs**
- `annotation_of_enh_for_potential_regulate_UpGene.bed`: Enahncers annotation of potentially regulated up-regulated genes
- `enh_expr_signal_for_UPgene.tab`: Tab-delimited signal table for enhancers
- `enh_sortedRegions_for_UPgene.bed`: Sorted enhancer regions
- `enh_UpGene_within_1mb.txt`: Enhancer-gene pairs for up-regulated genes within 1Mb
**Parameters**
- `Log2FC > 1`: Log2FC cutoffs for enhancer TPM signal

**Outputs**
- `SharedUPenhancer.txt`: List of shared upenhancers 
- `SharedUPgene_INenhancer.txt`: List of shared upregulated genes associated with shared upenhancers

**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Get_shared_Down_Enhancer_and_shared_Down_gene**

**Description**

Identifies enhancers that show consistent accessibility changes and are associated with downregulated genes across samples, creating lists of shared downenhancers and their associated genes.

**Inputs**
- `annotation_of_enh_for_potential_regulate_DownGene.bed`: Enahncers annotation of potentially regulated down-regulated genes
- `enh_expr_signal_for_DOWNgene.tab`: Tab-delimited signal table for enhancers
- `enh_sortedRegions_for_DOWNgene.bed`: Sorted enhancer regions
- `enh_DownGene_within_1mb.txt`: Enhancer-gene pairs for down-regulated genes within 1Mb
**Parameters**
- `Log2FC < -1`: Log2FC cutoffs for enhancer TPM signal

**Outputs**
- `SharedDOWNenhancer.txt`: List of shared downenhancers 
- `SharedDOWNgene_INenhancer.txt`: List of shared down-regulated genes associated with shared downenhancers

**Output directory**
- `./datasets/Enhancer_region_analysis/`

### **Enrichment_analysis_for_shared_gene**

**Description**

Performs Gene Set Enrichment Analysis (GSEA) for genes associated with shared enhancer regions, identifying biological processes and pathways that are significantly enriched in enhancer-associated gene sets.


**Inputs**
- `SharedUPgene_INenhancer.txt`
- `SharedDOWNgene_INenhancer.txt`

**Parameters**
**Parameters**
- Statistical significance threshold (p-value < 0.05 |NES| > 1)

**Outputs**
- `sharedUpGene_INenhancer_gsea_GOenrich.csv`: GO enrichment results for enhancer-associated upregulated genes
- `sharedUpGene_INenhancer_gsea_KEGGenrich.csv`: KEGG enrichment results for enhancer-associated upregulated genes
- `sharedDownGene_INenhancer_gsea_GOenrich.csv`: GO enrichment results for enhancer-associated downregulated genes
- `sharedDownGene_INenhancer_gsea_KEGGenrich.csv`: KEGG enrichment results for enhancer-associated downregulated genes
- `sharedUpGene_INenhancer_gsea_GOenrich.rds`: GO enrichment results for enhancer-associated upregulated genes(rds format)
- `sharedUpGene_INenhancer_gsea_KEGGenrich.rds`: KEGG enrichment results for enhancer-associated upregulated genes(rds format)
- `sharedDownGene_INenhancer_gsea_GOenrich.rds`: GO enrichment results for enhancer-associated downregulated genes(rds format)
- `sharedDownGene_INenhancer_gsea_KEGGenrich.rds`: KEGG enrichment results for enhancer-associated downregulated genes(rds format)

**Note:** The rds format file is for user custom drawing GSEAplot

**Output directory**
- `./datasets/Enhancer_region_analysis/Enrichment_analysis/Up`
- `./datasets/Enhancer_region_analysis/Enrichment_analysis/Down`
### **Extract_enhancer_sequences**

**Description**

Extracts DNA sequences from shared enhancer regions for motif enrichment analysis, providing sequence data for identification of transcription factor binding sites within active enhancers.

**Inputs**
- `SharedUPenhancer.txt`
- `SharedDOWNenhancer.txt`
- `EnhancerAtlasv2.0_eRNA.bed`

**Parameters**
- Enhancer sequence extraction parameters
- Sequence length specifications

**Outputs**
- `Regulation_analysis/Up/SharedUP_enhancer_sequences.fa`: FASTA sequences for shared upenhancers
- `Regulation_analysis/Down/SharedDOWN_enhancer_sequences.fa`: FASTA sequences for shared downenhancers

**Output directory**
- `./datasets/Enhancer_region_analysis/Regulation_analysis/`

### **Motif_enrich**

**Description**

Performs transcription factor motif enrichment analysis on enhancer sequences using the ame tool of MEME software to identify transcription factors that are significantly enriched in accessible enhancer regions associated with differentially expressed genes.

**Inputs**
- `Regulation_analysis/Up/SharedUP_enhancer_sequences.fa`: FASTA sequences for shared upenhancers
- `Regulation_analysis/Down/SharedDOWN_enhancer_sequences.fa`: FASTA sequences for shared downenhancers
- `./datasets/JASPAR2024_CORE/*.meme`:JASPAR2024 motif database in MEME format


**Outputs**
- `Regulation_analysis/Up/motif_enrich_result/analysis_complete.txt`: Completion flag for upregulated gene enhancer motif analysis
- `Regulation_analysis/Down/motif_enrich_result/analysis_complete.txt`: Completion flag for downregulated gene enhancer motif analysis

**Output directory**
- `./datasets/Enhancer_region_analysis/Regulation_analysis/`

### **Constructe_TF_Gene_network**

**Description**

Constructs transcription factor-gene regulatory networks based on enhancer motif enrichment results, linking transcription factors to genes through enhancer-mediated regulation.

**Inputs**
- `Regulation_analysis/Up/motif_enrich_result/`: Motif enrichment results
- `Regulation_analysis/Down/motif_enrich_result`:Motif enrichment results


**Outputs**
- `Regulation_analysis/Up/SharedUPenhancer_TF_gene_network.txt`: TF-gene network
- `Regulation_analysis/Down/SharedDOWNenhancer_TF_gene_network.txt`: TF-gene network

**Output directory**
- `./datasets/Enhancer_region_analysis/Regulation_analysis/`

### **Extract_enriched_motifs**

**Description**

Extracts and summarizes significantly enriched transcription factor motifs from enhancer motif enrichment analysis, providing lists of TFs potentially regulating genes through enhancer elements.

**Inputs**
- `Regulation_analysis/Up/motif_enrich_result/`: Motif enrichment results
- `Regulation_analysis/Down/motif_enrich_result`:Motif enrichment results

**Parameters**
- `FDR < 0.05`:FDR cutoffs for significance

**Outputs**
- `Regulation_analysis/Up/enriched_tfs_list.txt`: List of enriched TFs in upregulated gene enhancers
- `Regulation_analysis/Down/enriched_tfs_list.txt`: List of enriched TFs in downregulated gene enhancers

**Output directory**
- `./datasets/Enhancer_region_analysis/Regulation_analysis/`

### **Plot_TF_heatmap**

**Description**
Generates heatmap visualizations showing expression patterns of transcription factors enriched in enhancer regions, providing insights into TF activity in enhancer-mediated gene regulation.

**Inputs**
- `Regulation_analysis/Up/enriched_tfs_list.txt`: List of enriched TFs for upregulated genes
- `Regulation_analysis/Down/enriched_tfs_list.txt`: List of enriched TFs for downregulated genes
- `stringtiefile/experiment_vs_ctr_DEG.csv`:RNA-seq differential expression results

**Parameters**
- `Log2FC > 1`: Log2FC cutoffs for TF differential expression
- `FDR < 0.05`: FDR cutoffs for TF differential expression significance

**Outputs**
- `Regulation_analysis/Up/shared_UpEnhancer_tfs_heatmap_RNAseqExpr.pdf`: TF expression heatmap for enhancer-associated upregulated gene analysis
- `Regulation_analysis/Down/shared_DownEnhancer_tfs_heatmap_RNAseqExpr.pdf`: TF expression heatmap for enhancer-associated downregulated gene analysis

**Output directory**
- `./datasets/Enhancer_region_analysis/Regulation_analysis/`