required_packages <- c("Seurat", "Signac", "dplyr", "ggplot2", "GenomicRanges", "future", "patchwork", "hdf5r",
                       "readr", "pheatmap", "ggrepel", "LSD", "MASS", "ensembldb", "rtracklayer")
lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

counts_h3 <- Read10X_h5("..../filtered_feature_bc_matrix_h3.h5")
metadata_h3 <- read.csv(
  file = "..../per_barcode_metrics_h3.csv",
  header = TRUE,
  row.names = 1
)
fragpath_h3 <- "..../atac_fragments_h3.tsv.gz"

#Annotation
genome.h3 <- import("..../Danio_rerio.GRCz11.113.filtered.gtf")
genome(genome.h3) <- "GRCz11"

#create Seurat object
H3.data <- CreateSeuratObject(counts = counts_h3$`Gene Expression`, 
                              assay = 'RNA',
                              project = 'H3.3 KD',
                              meta.data = metadata_h3)

H3.data[["percent.mt"]] <- PercentageFeatureSet(H3.data, pattern = "^MT-")
atac_counts.h3 <- counts_h3$Peaks
grange.counts.h3 <- StringToGRanges(rownames(atac_counts.h3), sep = c(":", "-"))
grange.use.h3 <- seqnames(grange.counts.h3) %in% standardChromosomes(grange.counts.h3)
atac_counts.h3 <- atac_counts.h3[as.vector(grange.use.h3), ]

H3.data[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts.h3,
  sep = c(":", "-"),
  genome = genome(genome.h3),
  fragments = fragpath_h3,
  min.cells = 10,
  annotation = genome.h3
)

#QC
DefaultAssay(H3.data) <- "ATAC"
# compute nucleosome signal score per cell
H3.data <- NucleosomeSignal(object = H3.data)

# compute TSS enrichment score per cell
H3.data <- TSSEnrichment(object = H3.data, fast = FALSE)
DensityScatter(H3.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

#peak calling
h3.peaks <- system('bash -c ".../atac_fragments_h3.tsv.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n H33_KD --outdir ..../H3.3_sample/Peak_calling" ', 
                   wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)
peak.path_h3 <- "..../H33_KD_peaks.narrowPeak"
h3.peaks <- rtracklayer::import(peak.path_h3, format = "narrowPeak")
h3.peaks <- keepStandardChromosomes(h3.peaks, pruning.mode = 'coarse')

macs_count.h3 <- FeatureMatrix(fragments = Fragments(H3.data),
                               features = h3.peaks,
                               cells = colnames(H3.data))
H3.data[['peaks']] <- CreateChromatinAssay(
  counts = macs_count.h3,
  sep = c(":", "-"),
  genome = genome(genome.h3),
  fragments = fragpath_h3,
  min.cells = 10,
  annotation = genome.h3
)
DefaultAssay(H3.data) <- "peaks"
# compute nucleosome signal score per cell
H3.data <- NucleosomeSignal(object = H3.data)

# compute TSS enrichment score per cell
H3.data <- TSSEnrichment(object = H3.data, fast = FALSE)
DensityScatter(H3.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

saveRDS(H3.data, file = "...../sox10ncc_h33_kd.rds")

