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
library(BSgenome.Drerio.UCSC.danRer11)

set.seed(1234)

counts_control <- Read10X_h5("E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample/bin/filtered_feature_bc_matrix.h5")
metadata_control <- read.csv(
  file = "E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample/bin/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
fragpath_control <- "E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample/bin/atac_fragments.tsv.gz"

#Annotation
genome.control <- import("E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Reference_genome/Danio_rerio.GRCz11.113.filtered.gtf")
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
control.peaks <- system('bash -c "/home/ayush/.local/bin/macs2 callpeak -t /mnt/d/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample/Control_NewGRC/atac_fragments.tsv.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n Control --outdir /mnt/d/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample//Peak_calling" ', 
                        wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)
peak.path_c <- "E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/Control_sample/Peak_calling/Control_peaks.narrowPeak"
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

