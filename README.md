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
