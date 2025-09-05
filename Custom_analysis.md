# **RAGER table of contents**
1. [Quick start](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_RNAseq/RNAseq_analysis.md)
3. [Preprocess ATACseq data](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Preprocess_ATACseq/ATACseq_analysis.md) 
4. [Joint analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Joint_analysis/Joint_analysis.md)
5. [Custom analysis](https://github.com/yjliu15924/RAGER/blob/main/RAGER_github/Scripts/Custom_analysis/Custom_analysis.md)

## **Quick Start**

Welcome to RAGER! This section will guide through the essential steps to configure and run custom genes analysis pipeline.

Note that the custom analysis is based on the `preprocess RNAseq data`, `preprocess ATACseq data` and `Joint analysis`.
### **Step 1: Open the Configuration File**

Open the `config.yaml` file in the species-specific directory:

```bash
cd ~/PROJECT/RAGER  #Change the directory to PROJECT/RAGER as the current directory
#mouse
vim ./mouse/scripts/snakemake/Custom_genes_analysis/config.yaml  
#human
vim ./human/scripts/snakemake/Custom_genes_analysis/config.yaml 
# you can use any text editor 
```

### **Step 2: Modifying the Configuration File**

The `config.yaml` file contains all parameters needed to customize RAGER for custom genes analysis. Below is a detailed explanation of each section:

#### **Input Files Section**
Specify the input files for custom genes analysis:

```yaml
input_files:
  custom_genes: "datasets/Custom_genes_analysis/custom_genes.txt"    # Custom gene list file
  deg_csv: "datasets/RNAseq/stringtiefile/A_vs_N_DEG.csv"          # Differential expressed genes file
```

**How to modify**: 
- Update `custom_genes` path to point to your gene list file (one gene symbol per line)
- Update `deg_csv` path to match your RNA-seq differential expression analysis output
- Ensure both files are accessible and properly formatted

#### **Paths Section**
This section defines the directory structure for custom genes analysis:

```yaml
paths:
  dir: "datasets/Custom_genes_analysis"                    # Main analysis directory
  scripts_dir: "scripts/snakemake/Custom_genes_analysis/"  # Scripts directory
  
  jaspar_dir: "datasets/JASPAR2024_CORE"                  # JASPAR database directory
  tf_mapping_file: "datasets/Mus_musculus_TF.txt"         # Transcription factor mapping file
```

**How to modify**: 
- All directories will be created by the pipeline if they don't exist
- Update `tf_mapping_file` for different species (e.g.`datasets/Homo_sapiens_TF.txt` for human)

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
- These parameters will be used to filter the differential expression results for enrichment TFs.

### **Step 3: Prepare Custom Gene List**

Create your custom gene list file:

```bash
# Create the custom genes file
#mouse
mkdir -p ./mouse/datasets/Custom_genes_analysis
vim ./mouse/datasets/Custom_genes_analysis/custom_genes.txt

#human
mkdir -p ./human/datasets/Custom_genes_analysis
vim ./human/datasets/Custom_genes_analysis/custom_genes.txt
```

The custom genes file should contain one gene symbol per line:
```
AC016739
ACTR2
AHSP
ALAS2
APOE
ATP6V1E1
```

**File format requirements**:
- One gene symbol per line
- No headers or additional columns
- UTF-8 text encoding

### **Step 4: Validate Configuration**

Before running the pipeline, ensure:

1. Custom gene list file exists and is properly formatted
2. DEG CSV file from RNA-seq analysis is accessible
3. JASPAR database and TF mapping files are available
4. Differential expression thresholds are appropriate for your analysis
5. All file paths are correct and accessible
6. Output directory permissions allow write access

### **Step 5: Run the Pipeline**

Once configuration and custom gene list is ready, execute:

```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory
snakemake --snakefile ./scripts/snakemake/Custom_genes_analysis/Custom_gene_analysis.py --configfile ./scripts/snakemake/Custom_genes_analysis/config.yaml -j 10

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
snakemake --snakefile ./scripts/snakemake/Custom_genes_analysis/Custom_gene_analysis.py --configfile ./scripts/snakemake/Custom_genes_analysis/config.yaml -j 10
```

### **Expected Outputs**

The pipeline will generate analysis results in the `datasets/Custom_genes_analysis/` directory, including:

- Filtered differential expression results for your custom genes
- Transcription factor binding analysis
- Regulatory network analysis
- Visualization plots and summary statistics

**Note**: This analysis pipeline focuses on your specific genes of interest and provides targeted regulatory insights based on the custom gene list you provide.

## **List of processes**
- [Convert gene ID](#convert_customgene_to_ensemlid)
- [Enrichment_analysis](#)
- [Regulatory_analysis](#)
  - [Extract_custom_genes_promoterSeq](#)
  - [Motif_enrich](#)
  - [Constructe_TF_Gene_network](#)
  - [Extract_enriched_motifs](#)
  - [Plot_TF_heatmap](#)
- [Custom_plot](#custom_plot)

## **Preprocessing**
### **convert_customGene_to_ensemlID**

**Description**

Converts custom gene names/symbols to Ensembl gene IDs for standardized downstream analysis.

**Inputs**
- `custom_genes.txt`: List of custom gene names/symbols

**Outputs**
- `custom_genes_ensemblID.txt`: Custom genes converted to Ensembl IDs

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **GSEA_enrichment**

**Description**

Performs Gene Set Enrichment Analysis (GSEA) for custom genes using differential expression data, identifying enriched GO biological processes and KEGG pathways.

**Inputs**
- `custom_genes_ensemblID.txt`: Custom genes with Ensembl IDs
- `RNA_GDC_vs_RNA_ctr_DEG.csv`: Differential expression gene data

**Outputs**
- `custom_gene_gsea_GOenrich.csv`: GO enrichment analysis results
- `custom_gene_gsea_KEGGenrich.csv`: KEGG pathway enrichment analysis results
- `custom_gene_gsea_GOenrich.rds`: GO enrichment analysis results(rds format)
- `custom_gene_gsea_KEGGenrich.rds`: KEGG pathway enrichment analysis results(rds format)

**Note:** The rds format file is for user custom drawing GSEAplot

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **extr_promoterSeq**

**Description**

Extracts promoter sequences for custom genes in FASTA format for downstream motif analysis.

**Inputs**
- `custom_genes_ensemblID.txt`: Custom genes with Ensembl IDs


**Outputs**
- `promoter_seqs.fa`: Extracted promoter sequences in FASTA format

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **motif_enrich**

**Description**

Performs transcription factor binding motif enrichment analysis on promoter sequences using JASPAR motif database.

**Inputs**
- `promoter_seqs.fa`: Promoter sequences in FASTA format
- `JASPAR2024_CORE/`: JASPAR motif database directory

**Outputs**
- `motif_enrich_result`:Transcription factor binding motif enrichment result directory

**Output directory**
- `./datasets/Custom_genes_analysis/motif_enrich_result/`

---

### **TF_Gene_network**

**Description**

Constructs transcription factor-gene regulatory network based on motif enrichment analysis results.

**Inputs**
- `motif_enrich_result`: Transcription factor binding motif enrichment result directory
**Outputs**
- `TF_gene_network.txt`: Transcription factor-gene regulatory network

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **extract_enriched_motifs**

**Description**

Extracts significantly enriched transcription factor motifs and maps them to corresponding transcription factor names.

**Inputs**
- `motif_enrich_result`: Transcription factor binding motif enrichment result directory
- `Homo_sapiens_TF.txt`: Human transcription factor mapping file

**Parameters**
- `FDR < 0.05`:FDR cutoffs for significance

**Outputs**
- `enriched_tfs_list.txt`: List of significantly enriched transcription factors

**Output directory**
- `./datasets/Custom_genes_analysis/`

---

### **TF_heatmap**

**Description**

Generates expression heatmap for enriched transcription factors using differential expression data with specified fold change and significance thresholds.

**Inputs**
- `enriched_tfs_list.txt`: List of significantly enriched transcription factors
- `RNA_GDC_vs_RNA_ctr_DEG.csv`: Differential expression gene data

**Parameters**
- `log2fc_up`:  Log2 fold change threshold for upregulated genes
- `log2fc_down`: Log2 fold change threshold for downregulated genes
- `padj`: 0.05 (Adjusted p-value threshold for statistical significance)

**Outputs**
- `TF_expr_heatmap.pdf`: Expression heatmap of enriched transcription factors

**Output directory**
- `./datasets/Custom_genes_analysis/`

## **Custom_plot**
In this section, please change the directory to the corresponding species
```bash
#mouse
cd ~/PROJECT/RAGER/mouse  #Change the directory to PROJECT/RAGER/mouse as the current directory

#human
cd ~/PROJECT/RAGER/human  #Change the directory to PROJECT/RAGER/human as the current directory
```

### **GSEA_plot**

**Description**

Generates GSEA (Gene Set Enrichment Analysis) plots for specific GO terms or KEGG pathway from previously computed GSEA results, creating publication-ready visualizations.

**Inputs**
- `GO term`
  - `gseaGO_results.rds`: GSEA GO enrichment analysis results in RDS format
  - `GO_term_ID`: Specific GO term identifier (e.g., GO:0042127)
- `KEGG pathway`
  - `gseaKEGG_results.rds`: GSEA KEGG enrichment analysis results in RDS format
  - `KEGG_pathway_ID`: Specific KEGG pathway identifier (e.g., hsa04010)


**Outputs**
- `GSEA_plot_for_[GO_term].pdf`: GSEA enrichment plot for the specified GO term

**Command line**
```bash
Rscript --vanilla ./scripts/custom_plot/GSEA_plot.R [input_rds_file] [GO_term_ID] [output_pdf_file]

#Examples
Rscript --vanilla ./scripts/custom_plot/GSEA_plot.R ./datasets/Promoter_region_analysis/Enrichment_analysis/Up/sharedUpGene_InPromote_gsea_GOenrich.rds GO:0042127 ./GSEA_plot_for_GO_0042127.pdf
```

---

### **GSEA_network_to_cytoscape**

**Description**

Converts GSEA enrichment results and differential expression data into Cytoscape-compatible network format for pathway visualization and network analysis.

**Inputs**
- `Selected_gsea_enrich.csv`: Selected GSEA enrichment analysis results in CSV format
- `A_vs_N_DEG.csv`: Differential expression gene data

**Outputs**
- `GSEA_enrichment.graphml`: Cytoscape network file in GraphML format

**Command line**
```bash
Rscript --vanilla ./scripts/custom_plot/GSEA_network_to_cytoscape_demo.R [gsea_results.csv] [DEGgenes_data.csv] [output.graphml]

#Examples
Rscript --vanilla ./scripts/custom_plot/GSEA_network_to_cytoscape_demo.R ./datasets/Selected_gseaGO_enrich.csv ./Datasets/RNAseq/stringtiefile/A_vs_N_DEG.csv ./GOpathway_enrichment.graphml
```
**CRITICAL:**
1. `Selected_gsea_enrich.csv `:select the path you are interested in and need to keep the original format.
2. `Selected_gsea_enrich.csv` you should choose from one result, for example: sharedUpGene_InPromote_gsea_GOenrich.csv
3. `GSEA_enrichment.graphml`:This result you need to import Cytoscape software, where you can make the visualizations you want.For details, please refer to [Cytoscape](https://cytoscape.org/welcome.html).

### **Motif_enrich_to_cytoscape**

**Description**

Converts transcription factor-gene regulatory networks to Cytoscape-compatible format for network visualization.

**Inputs**
- `Shared_UPgene_INpromote_TF_gene_network.txt`: Transcription factor-gene regulatory network data
- `TF_list`: Comma-separated list of specific transcription factors of interest

**Outputs**
- `tf_gene_network.graphml`: Cytoscape network file in GraphML format


**Command line**
```bash
Rscript --vanilla ./scripts/custom_plot/Motif_enrich_tocytoscape_demo.R [tf_gene_network.txt] [output.graphml] [TF1,TF2,TF3,...]

#Examples
Rscript --vanilla ./scripts/custom_plot/Motif_enrich_tocytoscape_demo.R ./Datasets/Promoter_region_analysis/Regulation_analysis/Up/Shared_UPgene_INpromote_TF_gene_network.txt ./tf_gene_network.graphml KLF17,SNAI1,OTX2,KLF14,ZEB1
```
**CRITICAL:**
1. Before selecting TF, you need to look at ``enriched_tfs_list.txt`` and select TFs whose binding motifs are significantly enriched.
2. `TF_list`:select the Tfs you are interested in.
3. `tf_gene_network.graphml`:This result you need to import Cytoscape software, where you can make the visualizations you want.For details, please refer to [Cytoscape](https://cytoscape.org/welcome.html).

---