
# **RAGER table of contents**
1. [Quick start](https://github.com/bioinfo202408/RAGER/blob/main/README.md#quick-start)
2. [Preprocess RNAseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_RNAseq_data.md)
3. [Preprocess ATACseq data](https://github.com/bioinfo202408/RAGER/blob/main/Preprocess_ATACseq_data.md)
4. [Joint analysis](https://github.com/bioinfo202408/RAGER/blob/main/Joint_analysis.md)
5. [Custom analysis](https://github.com/bioinfo202408/RAGER/blob/main/Custom_analysis.md)
6. [UI](https://github.com/bioinfo202408/RAGER/blob/main/RAGER_UI.md)
7. [Unmapped reads module](https://github.com/bioinfo202408/RAGER/blob/main/Unmapped-reads_module.md)

# **Unmapped reads module**

This document describes the optional unmapped-read modules for RNA-seq and ATAC-seq in RAGER. These modules are designed to provide additional interpretation for reads that do not align to the reference genome, especially in species with incomplete or fragmented genome assemblies.

## **List of processes**
- [RNAseq_unmapped_reads_module](#rnaseq_unmapped_reads_module)
  - [Download reference files](#download-reference-files-for-rnaseq)
  - [Run the RNAseq unmapped module](#run-the-rnaseq-unmapped-module)
  - [RNAseq unmapped module outputs](#rnaseq-unmapped-module-outputs)
- [ATACseq_unmapped_reads_module](#atacseq_unmapped_reads_module)
  - [Download reference files](#download-reference-files-for-atacseq)
  - [Run the ATACseq unmapped module](#run-the-atacseq-unmapped-module)
  - [ATACseq unmapped module outputs](#atacseq-unmapped-module-outputs)

---

# RNAseq_unmapped_reads_module

## **Description**

This optional module extracts unmapped RNA-seq reads from the original alignment results and performs de novo transcript assembly followed by sequence annotation. It can be used to examine whether a portion of the unmapped RNA-seq reads may correspond to transcripts from incompletely assembled or incompletely annotated genomic regions.

## **Download reference files for RNAseq**

```bash
cd ~/PROJECT/RAGER/triticum_aestivum
mkdir -p ./reference/{cdna,protein}

# cdna
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz -P ./reference/cdna
gunzip ./reference/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz

# protein
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/triticum_aestivum/pep/Triticum_aestivum.IWGSC.pep.all.fa.gz -P ./reference/protein
gunzip ./reference/protein/Triticum_aestivum.IWGSC.pep.all.fa.gz
````

## **Run the RNAseq unmapped module**

### **Paired-end RNA-seq**

```bash
python3 ./scripts/unmapped-reads_module/rnaseq_unmapped_module.py \
  --input ./datasets/RNAseq/hisat2file/SRRxxxxxx/accepted_hits.sam \
  --sample SRRxxxxxx \
  --layout PE \
  --outdir ./datasets/RNAseq/unmapped_module \
  --threads 10 \
  --max-memory 100G \
  --min-contig-length 200 \
  --cdna-db ./reference/cdna/Triticum_aestivum.IWGSC.cdna.all.fa \
  --protein-db ./reference/protein/Triticum_aestivum.IWGSC.pep.all.fa
```

### **Single-end RNA-seq**

```bash
python3 ./scripts/unmapped-reads_module/rnaseq_unmapped_module.py \
  --input ./datasets/RNAseq/hisat2file/sample1/accepted_hits.sam \
  --sample sample1 \
  --layout SE \
  --outdir ./datasets/RNAseq/unmapped_module \
  --threads 10 \
  --max-memory 40G \
  --cdna-db ./reference/cdna/Triticum_aestivum.IWGSC.cdna.all.fa
```

## **RNAseq unmapped module outputs**

```bash
./datasets/RNAseq/unmapped_module/
├── SRRxxxxxx.unmapped_module.log
├── 01_extract_unmapped/
│   ├── SRRxxxxxx.unmapped.R1.fastq.gz
│   ├── SRRxxxxxx.unmapped.R2.fastq.gz
│   └── SRRxxxxxx.unmapped.singletons.fastq.gz
├── 02_trinity/
│   └── SRRxxxxxx_trinity/
│       └── Trinity.fasta
└── 03_annotation/
    └── SRRxxxxxx/
        ├── SRRxxxxxx.blastn.tsv
        ├── SRRxxxxxx.blastn.besthit.tsv
        ├── SRRxxxxxx.blastx.tsv
        └── SRRxxxxxx.blastx.besthit.tsv
```

### **SRRxxxxxx.unmapped_module.log**

**Description**

This log file records the commands executed by the RNAseq unmapped-read module and can be used to inspect runtime messages and possible errors.

### **01_extract_unmapped/**

**Description**

This directory contains unmapped RNA-seq reads extracted from the original alignment results.

**Outputs**

* `*.unmapped.R1.fastq.gz`: Unmapped read1 sequences from paired-end libraries
* `*.unmapped.R2.fastq.gz`: Unmapped read2 sequences from paired-end libraries
* `*.unmapped.singletons.fastq.gz`: Unpaired singleton reads generated during unmapped-read extraction

### **02_trinity/**

**Description**

This directory contains the de novo assembly results generated from unmapped RNA-seq reads.

**Outputs**

* `Trinity.fasta`: Assembled transcript contigs from unmapped RNA-seq reads

### **03_annotation/**

**Description**

This directory contains annotation results for assembled contigs based on similarity searches against cDNA and protein reference databases.

**Outputs**

* `*.blastn.tsv`: Full BLASTN results against the cDNA reference
* `*.blastn.besthit.tsv`: Best BLASTN hit retained for each contig
* `*.blastx.tsv`: Full BLASTX results against the protein reference
* `*.blastx.besthit.tsv`: Best BLASTX hit retained for each contig

---

# ATACseq_unmapped_reads_module

## **Description**

This optional module extracts unmapped ATAC-seq reads from the original alignment results and classifies them sequentially against chloroplast, mitochondrial, and repeat/transposable element reference sequences. Reads that remain unmatched after these steps are retained as an unclassified fraction. This module provides a simple interpretation of unmapped ATAC-seq reads.

## **Download reference files for ATACseq**

```bash
cd ~/PROJECT/RAGER/triticum_aestivum
mkdir -p ./reference/{organelle,repeat_TE}

# chloroplast
wget -c "https://www.ncbi.nlm.nih.gov/nuccore/NC_002762.1?report=fasta" --no-check-certificate -O ./reference/organelle/chloroplast.fa

# mitochondrion
wget -c "https://www.ncbi.nlm.nih.gov/nuccore/NC_036024.1?report=fasta" --no-check-certificate -O ./reference/organelle/mitochondria.fa
# If the download fails, please download the FASTA files manually from NCBI.

# repeat_TE
wget -c https://trep-db.uzh.ch/downloads/trep-db_complete_Rel-25.fasta.gz -P ./reference/repeat_TE
gunzip ./reference/repeat_TE/trep-db_complete_Rel-25.fasta.gz
```

## **Run the ATACseq unmapped module**

```bash
python3 ./scripts/unmapped-reads_module/atacseq_unmapped_module.py \
  --input ./datasets/ATACseq/bowtie2file/SRRxxxxxx/accepted_hits.sam \
  --sample SRRxxxxxx \
  --outdir ./datasets/ATACseq/unmapped_module \
  --threads 10 \
  --chloroplast-ref ./reference/organelle/chloroplast.fa \
  --mitochondria-ref ./reference/organelle/mitochondria.fa \
  --te-ref ./reference/repeat_TE/trep-db_complete_Rel-25.fasta
```

## **ATACseq unmapped module outputs**

```bash
./datasets/ATACseq/unmapped_module/
├── 01_extract_unmapped/
│   └── SRRxxxxxx.unmapped.fa
├── 02_classification/
│   ├── chloroplast/
│   ├── mitochondrial/
│   └── repeat_TE/
├── 03_summary/
│   ├── SRRxxxxxx.unmapped_classification.long.tsv
│   └── SRRxxxxxx.unmapped_classification.wide.tsv
└── SRRxxxxxx.atac_unmapped_module.log
```

### **01_extract_unmapped/**

**Description**

This directory contains unmapped ATAC-seq reads extracted from the original alignment results.

**Outputs**

* `*.unmapped.fa`: FASTA file containing unmapped ATAC-seq reads

### **02_classification/**

**Description**

This directory contains the intermediate classification results from sequential remapping of unmapped ATAC-seq reads.

#### **chloroplast/**

**Description**

This directory contains alignment results against the chloroplast genome reference.

#### **mitochondrial/**

**Description**

This directory contains alignment results against the mitochondrial genome reference.

#### **repeat_TE/**

**Description**

This directory contains alignment results against the repeat/transposable element reference.

### **03_summary/**

**Description**

This directory contains the final summary tables for the unmapped ATAC-seq read classification results.

**Outputs**

* `*.unmapped_classification.long.tsv`: Long-format summary table containing read counts and percentages for each category
* `*.unmapped_classification.wide.tsv`: Wide-format summary table containing the total number of unmapped reads and the counts and percentages of chloroplast, mitochondrial, repeat/TE, and unclassified reads

### **SRRxxxxxx.atac_unmapped_module.log**

**Description**

This log file records the commands executed by the ATACseq unmapped-read module and can be used to inspect runtime messages and possible errors.

```
