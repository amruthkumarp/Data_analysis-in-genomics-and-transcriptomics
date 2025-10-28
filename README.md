# README — Multi-Omics Data Analysis Codes
---
**Purpose:**  
This repository contains only the **code** for multiple omics pipelines.  
  

---

## 1️⃣ Microarray Analysis — `oligo` package


**Techniques included:**
- Reading and processing `.CEL` files  
- Background correction and RMA normalization  
- Probe summarization  
- Differential expression analysis using `limma`  
- Annotation mapping to gene symbols  


## 2️⃣ Bulk RNA-seq — `DESeq2` and `edgeR`

 
**Techniques included:**
- Count matrix import (featureCounts)  
- Normalization (size factors, CPM, TMM)  
- Differential expression via `DESeq2` and `edgeR`  
- Volcano and MA plots for visualization  



---

## 3️⃣ ChIP-seq — `BWA` + `MACS2` peak calling


**Techniques included:**
- Alignment using `bwa mem`  
- Sorting and indexing BAM files with `samtools`  
- Peak calling using `macs2 callpeak`  
- Optional replicate merging and differential peak detection  


---

## 4️⃣ Alternative Splicing Analysis — `scseq2'
  
**Techniques included:**
- Junction count quantification from aligned BAM files  
- Identification of skipped exons, A5SS, A3SS, RI, MXE events  
- Differential splicing analysis between conditions  

---


## 6️⃣ Bisulfite Sequencing — `Bismark`


**Techniques included:**
- Alignment of bisulfite-treated reads to reference genome  
- Deduplication and methylation extraction  
- Methylation percentage calculation and report generation  


---

## 7️⃣ Variant Calling — `GATK` pipeline
  
**Techniques included:**
- Alignment (Bowtie2)  
- Sorting, duplicate marking, base recalibration  
- Variant calling using `GATK HaplotypeCaller`  
- Joint genotyping and VCF output  
---



**Maintainer:** *Amruthkumar P*  
**Last updated:** 2025-10-28  
**Scope:** Code-only version for microarray, RNA-seq, ChIP-seq, splicing, scRNA-seq, bisulfite, and variant calling pipelines.  
