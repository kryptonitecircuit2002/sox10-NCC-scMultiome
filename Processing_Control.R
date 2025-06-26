library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(future)
library(patchwork)
library(hdf5r)
library(readr)
library(pheatmap)
library(ggrepel)
library(LSD)
library(MASS)
library(ensembldb)
setwd("...")
set.seed(1234)

counts_control <- Read10X_h5("..../filtered_feature_bc_matrix.h5")
metadata_control <- read.csv(
  file = ".../per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
fragpath_control <- "..../atac_fragments.tsv.gz"

#Annotation
genome.control <- import("D:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Reference_genome/Danio_rerio.GRCz11.113.filtered.gtf")
genome(genome.control) <- "GRCz11"

#create Seurat object
Control.data <- CreateSeuratObject(counts = counts_control$`Gene Expression`, 
                                   assay = 'RNA',
                                   project = 'Control',
                                   meta.data = metadata_control)

Control.data[["percent.mt"]] <- PercentageFeatureSet(Control.data, pattern = "^MT-")
atac_counts.c <- counts_control$Peaks
grange.counts.c <- StringToGRanges(rownames(atac_counts.c), sep = c(":", "-"))
grange.use.c <- seqnames(grange.counts.c) %in% standardChromosomes(grange.counts.c)
atac_counts.c <- atac_counts.c[as.vector(grange.use.c), ]

Control.data[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts.c,
  sep = c(":", "-"),
  genome = genome(genome.control),
  fragments = fragpath_control,
  min.cells = 10,
  annotation = genome.control
)

#QC
DefaultAssay(Control.data) <- "ATAC"
# compute nucleosome signal score per cell
Control.data <- NucleosomeSignal(object = Control.data)

# compute TSS enrichment score per cell
Control.data <- TSSEnrichment(object = Control.data, fast = FALSE)
DensityScatter(Control.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

#Peak-calling
#peak calling
control.peaks <- system('bash -c ".../macs2 callpeak -t ../atac_fragments.tsv.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n Control --outdir .../Peak_calling" ', 
                        wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)
peak.path_c <- ".../Control_peaks.narrowPeak"
control.peaks <- rtracklayer::import(peak.path_c, format = "narrowPeak")
control.peaks <- keepStandardChromosomes(control.peaks, pruning.mode = 'coarse')

macs_count.control <- FeatureMatrix(fragments = Fragments(Control.data),
                                    features = control.peaks,
                                    cells = colnames(Control.data))
Control.data[['peaks']] <- CreateChromatinAssay(
  counts = macs_count.control,
  sep = c(":", "-"),
  genome = genome(genome.control),
  fragments = fragpath_control,
  min.cells = 10,
  annotation = genome.control
)
DefaultAssay(Control.data) <- "peaks"
Control.data <- NucleosomeSignal(object = Control.data)
Control.data <- TSSEnrichment(object = Control.data, fast = FALSE)
DensityScatter(Control.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

DensityScatter(Control.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = Control.data,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0) 

#Filter out low quality cells
Control.data  <- subset(
  x = Control.data,
  subset = nCount_ATAC < 150000 &
    nCount_RNA < 40000 &
    nucleosome_signal < 1.75 &
    TSS.enrichment < 10
)

#RNA ASsay normalization
Control.data <-  SCTransform(Control.data, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

#ATAC normalisation
DefaultAssay(Control.data) <- "peaks"
Control.data <- RunTFIDF(Control.data)
Control.data <- FindTopFeatures(Control.data, min.cutoff = "q5")
Control.data <- RunSVD(Control.data)
Control.data <- RegionStats(Control.data, genome = BSgenome.Drerio.UCSC.danRer11)
Control.data <- LinkPeaks(Control.data,
                          peak.assay = "peaks",
                          expression.assay = "SCT",
                          gene.id = TRUE) 
##ATAC assay
DefaultAssay(Control.data) <- "ATAC"
Control.data <- RunTFIDF(Control.data)
Control.data <- FindTopFeatures(Control.data, min.cutoff = "q10")
Control.data <- RunSVD(Control.data)
DensityScatter(Control.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) + ggtitle("TSS Fragment Enrichment in Control")

# PCA, clustering, and UMAP visualization
DefaultAssay(Control.data) <- "SCT"
Control.data <- RunPCA(Control.data, verbose = FALSE, reduction.name = "pca.control")
Control.data <- FindNeighbors(Control.data, reduction = "pca.control", dims = 1:30)
Control.data <- FindClusters(Control.data, resolution = 0.7)
Control.data <- RunUMAP(Control.data, reduction = "pca.control", reduction.name = "umap.rna.control", dims = 1:30, verbose = FALSE, spread = 0.32, min.dist = 0.35)
DimPlot(Control.data, reduction = "umap.rna.control", label = TRUE)
