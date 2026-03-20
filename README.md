# **RAGER: A User-Friendly Computational Platform for Integrated Analysis of RNA-Seq and ATAC-Seq Data**

## Introduction to RAGER:
**RAGER** is a computational platform, that integrates the popular bioinformatics tools in an automated thread for joint mining of RNA-seq and ATAC-seq data. RAGER facilitates integrative analysis of transcriptome and chromatin accessibility by providing an automated workflow that minimizes the need for bioinformatics expertise and significantly reduces processing time. We demonstrate RAGER's utility for novel biological discovery by characterizing the transcriptome and chromatin accessibility of two recently published datasets.

## Table of Contents
1. [Quick start](#quick-start)
2. [Preprocess RNAseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_RNAseq_data.md)
3. [Preprocess ATACseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_ATACseq_data.md) 
4. [Joint analysis](https://github.com/bioinfo202408/RAGER/blob/main/Joint_analysis.md)
5. [Custom analysis](https://github.com/bioinfo202408/RAGER/blob/main/Custom_analysis.md)
6. [UI](https://github.com/bioinfo202408/RAGER/blob/main/RAGER_UI.md)
7. [Unmapped reads module](https://github.com/bioinfo202408/RAGER/blob/main/Unmapped-reads_module.md)

# **Quick start**
## System requirements:
Some of the tools that RAGER uses, e.g. Hisat2 and Bowtie2 are very memory intensive programs. Therefore we recommend the following system requirements for RAGER:

### Minimal system requirements:
- **Operating System:** Linux (tested on CentOS Stream)  
- **Memory:** ≥ 16 GB RAM  
- **CPU:** ≥ 10 cores 
- **Storage:** ≥ 500 GB available space  

### Recommended system requirements:
We recommend that you have at least 1T of ram and at least a 20-core CPU if you want to run RAGER in multi-threaded mode (which will speedup the workflow significantly). Our own servers equipped with 192 CPU cores, 1.5 TB of RAM, 80 TB of storage, and running CentOS Stream.


## Getting Started - Install the Conda virtual environment and the necessary data: 

### Install virtual environment.
1. Install conda. Download the conda installer from https://www.anaconda.com/. 
```bash
bash <conda-installer-name>-latest-Linux-x86_64.sh
conda --version
```
(onda-installer-name)This is a placeholder. During actual use, it needs to be replaced with the specific installation package name.The conda version should appear if it has been installed correctly.

2. install git.
```bash
git --version
```
If it does not return the installed version of git, follow https://github.com/git-guides/install-git to install git.

3. Clone the branch from the GitHub repository which contains the RAGER code and the necessary files for this pipeline.
```bash
cd ~
mkdir PROJECT #you can name your PROJECT folder anything
cd PROJECT
git clone https://github.com/bioinfo202408/RAGER
```


4. Create RAGER virtual environment.
```bash
cd RAGER
conda env create -f rager.yml
```
We have provided a complete Conda environment file (**`rager.yml`**) in the GitHub repository, which includes all the dependencies required for RAGER. Users can download this file, optionally modify the environment name parameter (default `"rager"`), and create the environment using Conda. In addition, **MEME Suite AME** (v5.5.7) and the R package **`genekitr`** need to be installed separately.

 **1). MEME Suite AME**
 
**MEME Suite AME** is one of the tools required for analysis, and its installation requires root permission. Please refer to the following link for the installation guide: [MEME Suite AME (v5.5.7)](https://meme-suite.org/meme/meme_5.5.7/). Choose the appropriate installation method based on your environment.

---

 **2). Installing the R Package `genekitr`**
 
The R package **`genekitr`** is a crucial dependency for the analysis but cannot be installed via Conda. After setting up and activating the Conda environment, you need to install `genekitr` manually in R. Follow the steps below:

 **Step 1: Activate the Conda environment**
```bash
conda activate rager
```

 **Step 2: Open R or RStudio**
If you are working in a terminal, you can start R by simply typing:

```bash
R
```

This will enter the R interactive environment in the terminal interface.

Alternatively, if you are using RStudio, you can launch it from your graphical user interface (GUI) and open a new R script or console session.


 **Step 3): Install `genekitr` from CRAN**
Install the stable version of `genekitr` using the following command:

```r
install.packages("genekitr")
```

 **Step 4): Verify the installation**
After installation, ensure the package has been successfully installed:

```r
library(genekitr)
packageVersion("genekitr")
```


 **Step 5): Exit the R environment**
To exit R from the terminal or RStudio, run the following command:

```r
q()
```

R will prompt you to save the workspace:

- Type `"y"` and hit `Enter` to save the workspace.

---

## Download static reference files
RAGER is dependent on reference files which can be found for the supported species listed below: download link [hg38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/) [mm10](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/) [TAIR10](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/) [IWGSC](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/)


1. If you are analyzing human data,navigate to the `PROJECT/RAGER/human` folder.
```bash
cd ~/PROJECT/RAGER/human
#Download and unzip reference genome
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz -P ./reference

gunzip ./reference/GRCh38.p14.genome.fa.gz
#Index reference genome
bowtie2-build --threads 5 ./reference/GRCh38.p14.genome.fa ./reference/bowtie2index/GRCh38

hisat2-build -p 5 ./reference/GRCh38.p14.genome.fa ./reference/hisat2index/GRCh38
#Download and unzip annotation file
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz -P ./reference/geneanno

gunzip ./reference/geneanno/gencode.v44.annotation.gtf.gz
```
2. If you are analyzing mouse data,navigate to the `PROJECT/RAGER/mouse` folder.

```bash
cd ~/PROJECT/RAGER/mouse
#Download and unzip reference genome
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz -P ./reference/

gunzip ./reference/GRCm38.p6.genome.fa.gz

#Index reference genome
bowtie2-build --threads 5 ./reference/GRCm38.p6.genome.fa ./reference/bowtie2index/GRCm38

hisat2-build -p 5 ./reference/GRCm38.p6.genome.fa ./reference/hisat2index/GRCm38
#Download and unzip annotation file
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz -P ./reference/geneanno

gunzip ./reference/geneanno/gencode.vM25.annotation.gtf.gz
```

3. If you are analyzing arabidopsis thaliana data,navigate to the `PROJECT/RAGER/arabidopsis_thaliana` folder.

```bash
cd ~/PROJECT/RAGER/arabidopsis_thaliana
#Download and unzip reference genome
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -P ./reference

gunzip ./reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

#Index reference genome
bowtie2-build --threads 5 ./reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ./reference/bowtie2index/Arabidopsis_thaliana_TAIR10

hisat2-build -p 5 ./reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ./reference/hisat2index/Arabidopsis_thaliana_TAIR10
#Download and unzip annotation file
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gtf.gz -P ./reference/geneanno

gunzip ./reference/geneanno/Arabidopsis_thaliana.TAIR10.62.gtf.gz
```

4. If you are analyzing triticum aestivum data,navigate to the `PROJECT/RAGER/triticum_aestivum` folder.

```bash
cd ~/PROJECT/RAGER/triticum_aestivum
#Download and unzip reference genome
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz -P ./reference

gunzip ./reference/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

#Index reference genome
bowtie2-build --threads 5 ./reference/Triticum_aestivum.IWGSC.dna.toplevel.fa ./reference/bowtie2index/Triticum_aestivum_IWGSC

hisat2-build -p 5 ./reference/Triticum_aestivum.IWGSC.dna.toplevel.fa ./reference/hisat2index/Triticum_aestivum_IWGSC
#Download and unzip annotation file
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.62.gtf.gz -P ./reference/geneanno

gunzip ./reference/geneanno/Triticum_aestivum.IWGSC.62.gtf.gz

#Build R annotation packages (TxDb, OrgDb, and BSgenome)
cd reference
# ----------------------------------------------------------------------
# Build TxDb (SQLite) and OrgDb (gene annotation) packages
# ----------------------------------------------------------------------
Rscript build_TxDb_OrgDb.R ./geneanno/Triticum_aestivum.IWGSC.62.gtf ./

# Install the generated OrgDb package
R CMD build org.Taestivum.eg.db
R CMD INSTALL org.Taestivum.eg.db_1.0.tar.gz

# ----------------------------------------------------------------------
# Build BSgenome package (genome sequence)
# ----------------------------------------------------------------------
# Download the UCSC faToTwoBit utility (Linux 64‑bit)
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod 755 faToTwoBit

# Convert the FASTA genome to 2bit format (required by BSgenome)
./faToTwoBit Triticum_aestivum.IWGSC.dna.toplevel.fa Triticum_aestivum.IWGSC.dna.toplevel.2bit

# Run the BSgenome build script
Rscript build_BSgenome.R ./Triticum_aestivum.IWGSC.dna.toplevel.2bit ./

# Install the generated BSgenome package
R CMD build BSgenome.Taestivum.IWGSC.v1
R CMD INSTALL BSgenome.Taestivum.IWGSC.v1_1.0.0.tar.gz

# After successful installation, you can load the packages in an R session:
R
library(org.Taestivum.eg.db)
library(BSgenome.Taestivum.IWGSC.v1)
```

5. If you are analyzing your own reference genome (custom gene fragments), navigate to the corresponding species folder and build indexes.

```bash
#1. Navigate to the corresponding folder (choose ONE based on your species)
#Human
cd ~/PROJECT/RAGER/human
#Mouse
cd ~/PROJECT/RAGER/mouse
#Arabidopsis thaliana
cd ~/PROJECT/RAGER/arabidopsis_thaliana
#Triticum aestivum
cd ~/PROJECT/RAGER/triticum_aestivum

#2. Upload/put your reference genome fasta file into ./reference
#   Example file name: custom_reference.fa
#   (Make sure the file is in FASTA format and ends with .fa)

#3. Build indexes for your custom reference genome
#   Replace ./reference/custom_reference.fa with your real fasta file name
bowtie2-build --threads 5 ./reference/custom_reference.fa ./reference/bowtie2index/custom_reference

hisat2-build -p 5 ./reference/custom_reference.fa ./reference/hisat2index/custom_reference
```

## Download or generate RNA-seq and ATACseq datasets
RAGER has been verified to have excellent functionality in datasets GSE85632 and GSE261119.
Raw sequencing data were downloaded from the GEO database. From GSE85632 dataset, we selected specific SRR files (SRR4032350, SRR4032351, SRR4032352, SRR4032353,SRR4032269, SRR4032270, SRR4032271, SRR4032272) relevant to our analysis, resulting in a total of 16 fastq files.

From GSE261119 dataset, we selected specific SRR files (RNAseq:SRR28264346,SRR28264347,SRR28264348,SRR28264349;ATACseq:SRR28263042,SRR28263043,SRR28263044,SRR28263045) relevant to our analysis, resulting in a total of 16 fastq files.

SRR4032350-SRR4032353 is RNAseq data, please save fastqfile to `./datasets/RNAseq/fastqfile` directory. Use the SRR4032350 file as an example. 

This is a mouse data,navigate to the `PROJECT/RAGER/mouse` folder.
```bash
cd ~/PROJECT/RAGER/mouse #If you want to analyz human data,navigate to the ~/PROJECT/RAGER/human folder.
prefetch SRR4032350 -O ./datasets/rawdata

fastq-dump --split-files -O ./datasets/RNAseq/fastqfile ./datasets/rawdata/SRR4032350/SRR4032350.sra 
```


SRR4032269-SRR4032272 is ATACseq data, please save fastqfile to `./datasets/ATACseq/fastqfile` directory.Use the SRR4032269 file as an example.

```bash
prefetch SRR4032269 -O ./datasets/rawdata

fastq-dump --split-files -O ./datasets/ATACseq/fastqfile ./datasets/rawdata/SRR4032269/SRR4032269.sra 
```

**Users who have custom sequencing data only need to unzip the original data and rename it to a unified filename.**

**Note:** Place different sequencing data in their respective folders.

```bash
#mouse
#RNAseq data
gunzip -c test_R1.fq.gz > ~/PROJECT/RAGER/mouse/datasets/RNAseq/fastqfile/test_1.fastq
gunzip -c test_R2.fq.gz > ~/PROJECT/RAGER/mouse/datasets/RNAseq/fastqfile/test_2.fastq

#ATACseq data
gunzip -c test2_R1.fq.gz > ~/PROJECT/RAGER/mouse/datasets/ATACseq/fastqfile/test2_1.fastq
gunzip -c test2_R2.fq.gz > ~/PROJECT/RAGER/mouse/datasets/ATACseq/fastqfile/test2_2.fastq

#human
#RNAseq data
gunzip -c test_R1.fq.gz > ~/PROJECT/RAGER/human/datasets/RNAseq/fastqfile/test_1.fastq
gunzip -c test_R2.fq.gz > ~/PROJECT/RAGER/human/datasets/RNAseq/fastqfile/test_2.fastq

#ATACseq data
gunzip -c test2_R1.fq.gz > ~/PROJECT/RAGER/human/datasets/ATACseq/fastqfile/test2_1.fastq
gunzip -c test2_R2.fq.gz > ~/PROJECT/RAGER/human/datasets/ATACseq/fastqfile/test2_2.fastq
```
In short, please unzip the raw data to *_1.fastq *_2.fastq format and place in their respective folders.
## **Comparative Features of Bioinformatics Tools for Epigenomics Data Analysis** 
| Classification          | Trait                                   | DiffBind | PRADA | TRAPLINE | esATAC | CoBRA | RAGER |
|-------------------------|-----------------------------------------|----------|-------|----------|--------|-------|-------|
| **Basic Information**   | Year                                    | 2011     | 2014  | 2016     | 2018   | 2021  | 2025  |
|                         | Workflow Engine                         | -        | -     | Galaxy   | -      | Snakemake | Snakemake |
|                         | Containerization Support                | -        | -     | -        | -      | ✓     | ✓     |
| **Data Preprocessing Module** | Standard Quality Control and Comparison | ✓        | ✓     | ✓        | ✓      | -     | ✓     |
|                         | Peak Calling                            | ✓        | -     | -        | ✓      | ✓     | ✓     |
| **Advanced Module**     | Joint Analysis                          | -        | -     | -        | -      | ✓     | ✓     |
|                         | Custom Analysis                         | -        | -     | -        | ✓      | ✓     | ✓     |
|                         | Motif Enrichment                        | ✓        | -     | -        | ✓      | ✓     | ✓     |
|                         | Differential Analysis                   | ✓        | ✓     | ✓        | -      | -     | ✓     |
|                         | Chromatin State Enrichment              | ✓        | -     | -        | ✓      | ✓     | ✓     |
| **User Experience**     | One-click Report Generation             | -        | ✓     | ✓        | -      | -     | ✓     |
|                         | Detailed instruction with case studies  | -        | ✓     | ✓        | -      | ✓     | ✓     |

---

### **The subsequent analysis steps should be referred to the readme files of each process.**














