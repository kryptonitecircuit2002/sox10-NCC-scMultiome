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

counts <- Read10X_h5(".../bin/filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = ".../bin/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
fragpath <- ".../bin/atac_fragments.tsv.gz"

#Annotation
genome <- import(".../Danio_rerio.GRCz11.113.filtered.gtf")
genome(genome) <- "GRCz11"


#create Seurat object
data <- CreateSeuratObject(counts = counts$`Gene Expression`, 
                                   assay = 'RNA',
                                   project = 'data',
                                   meta.data = metadata)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
atac_counts.c <- counts$Peaks
grange.counts.c <- StringToGRanges(rownames(atac_counts.c), sep = c(":", "-"))
grange.use.c <- seqnames(grange.counts.c) %in% standardChromosomes(grange.counts.c)
atac_counts.c <- atac_counts.c[as.vector(grange.use.c), ]

data[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts.c,
  sep = c(":", "-"),
  genome = genome(genome),
  fragments = fragpath,
  min.cells = 10,
  annotation = genome
)

#QC
DefaultAssay(data) <- "ATAC"
# compute nucleosome signal score per cell
data <- NucleosomeSignal(object = data)

# compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)
DensityScatter(data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

#Peak-calling
system('bash -c "...bin/macs2 callpeak -t /mnt/..../atac_fragments.tsv.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n data --outdir .../Peak_calling" ', 
                        wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)
peak.path_c <- ".../Peak_calling/data_peaks.narrowPeak"
peaks <- rtracklayer::import(peak.path_c, format = "narrowPeak")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')

macs_count <- FeatureMatrix(fragments = Fragments(data),
                                    features = peaks,
                                    cells = colnames(data))
data[['peaks']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = genome(genome),
  fragments = fragpath,
  min.cells = 10,
  annotation = genome
)
DefaultAssay(data) <- "peaks"
data <- NucleosomeSignal(object = data)
data <- TSSEnrichment(object = data, fast = FALSE)
DensityScatter(data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

saveRDS(data, file = "..../preprocessed_data.rds")

