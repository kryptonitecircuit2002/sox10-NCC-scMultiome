# Multiome Analysis Pipeline for Zebrafish Neural Crest Cells

This repository contains a computational workflow for the analysis of multimodal single-cell multiome datasets generated from sox10⁺ neural crest cells (NCCs). The pipeline uses Seurat v5.1.0 and Signac v1.15.0 for quality control, integration, dimensionality reduction, differential analysis and visualization of paired single-cell transcriptomic and chromatin accessibility data across experimental conditions.

## Experimental System

Samples were generated from Tg(sox10:GFP) zebrafish (Danio rerio) embryos at 24 hours post-fertilization (hpf).

## Sequencing Platform

- 10x Genomics Single Cell ATAC + Gene Expression (3' Gene Expression) Multiome Kit
- Illumina NextSeq 1000

## Input Files

The pipeline requires the following outputs generated using Cell Ranger ARC v8:

- `.h5` filtered feature-barcode matrix
- `.tsv.gz` fragments file
- `barcode_metrics.csv`

## Software

- R
- Seurat v5.1.0
- Signac v1.15.0

## Repository Structure

```text
├── Control
│   ├── 01_PreProcessing.R
│   ├── 02_Clustering.R
│   ├── 03_Links.R
│   ├── 04_Celltype Annotation.R
│   ├── 05_TF Motif Enrichment.R
│   ├── 06_Module Scoring.R
│   └── Pseudotime_TSCAN.R
│
├── H2AZ_knockdown
│   ├── 01_Preprocessing.R
│   ├── 02_Integration.R
│   ├── 03_Annotation.R
│   ├── 04_Module Scoring.R
│   └── 05_CytoTRACE.R
│
└── H33_knockdown
    ├── 01_Preprocessing.R
    ├── 02_Integration.R
    ├── 03_ATAC Integration.R
    ├── 04_Annotation.R
    ├── 05_Module Scoring.R
    ├── 06_CytoTRACE.R
    └── 07_TF_Motif_enrichment.R
```

## Analysis Workflow

### Control Dataset

- **01_PreProcessing.R** – Quality control, filtering, normalization, and multimodal object generation.
- **02_Clustering.R** – Dimensionality reduction, neighborhood graph construction, and clustering.
- **03_Links.R** – Peak-to-gene linkage analysis.
- **04_Celltype Annotation.R** – Cell type annotation using marker genes and accessibility profiles.
- **05_TF Motif Enrichment.R** – Motif enrichment analysis using chromatin accessibility data.
- **06_Module Scoring.R** – Gene regulatory module scoring.
- **Pseudotime_TSCAN.R** – Developmental trajectory inference using TSCAN.

### H2A.Z Knockdown

- **01_Preprocessing.R** – Data preprocessing and QC.
- **02_Integration.R** – Dataset integration.
- **03_Annotation.R** – Cell type annotation.
- **04_Module Scoring.R** – Module scoring analysis.
- **05_CytoTRACE.R** – Developmental potential analysis using CytoTRACE.

### H3.3 Knockdown

- **01_Preprocessing.R** – Data preprocessing and QC.
- **02_Integration.R** – Dataset integration.
- **03_ATAC Integration.R** – Integration of transcriptomic and chromatin accessibility modalities.
- **04_Annotation.R** – Cell type annotation.
- **05_Module Scoring.R** – Module scoring analysis.
- **06_CytoTRACE.R** – Developmental potential analysis using CytoTRACE.
- **07_TF_Motif_enrichment.R** – Transcription factor motif enrichment analysis.
