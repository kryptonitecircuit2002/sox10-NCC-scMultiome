required_packages <- c("Seurat", "Signac", "dplyr", "ggplot2", "GenomicRanges", "future", "patchwork", "hdf5r",
                       "readr", "pheatmap", "ggrepel", "LSD", "MASS", "ensembldb", "rtracklayer")
lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

counts_h2 <- Read10X_h5("...../filtered_feature_bc_matrix_h2.h5")

metadata_h2 <- read.csv(
  file = "..../per_barcode_metrics_h2.csv",
  header = TRUE,
  row.names = 1
)
fragpath_h2 <- "..../atac_fragments.tsv_h2.gz"

#Annotation
genome.h2 <- import("..../Danio_rerio.GRCz11.113.filtered.gtf")
genome(genome.h2) <- "GRCz11"
seqlevelsStyle(genome.h2) <- "UCSC"  

#create Seurat object
H2.data <- CreateSeuratObject(counts = counts_h2$`Gene Expression`, 
                              assay = 'RNA',
                              project = 'H2A.Z KD',
                              meta.data = metadata_h2)

H2.data[["percent.mt"]] <- PercentageFeatureSet(H2.data, pattern = "^MT-")
atac_counts.h2 <- counts_h2$Peaks
grange.counts.h2 <- StringToGRanges(rownames(atac_counts.h2), sep = c(":", "-"))
grange.use.h2 <- seqnames(grange.counts.h2) %in% standardChromosomes(grange.counts.h2)
atac_counts.h2 <- atac_counts.h2[as.vector(grange.use.h2), ]

H2.data[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts.h2,
  sep = c(":", "-"),
  genome = genome(genome.h2),
  fragments = fragpath_h2,
  min.cells = 10,
  annotation = genome.h2
)

#QC
DefaultAssay(H2.data) <- "ATAC"
# compute nucleosome signal score per cell
H2.data <- NucleosomeSignal(object = H2.data)

# compute TSS enrichment score per cell
H2.data <- TSSEnrichment(object = H2.data, fast = FALSE)
DensityScatter(H2.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

#peak calling
h2.peaks <- system('bash -c "..../atac_fragments.tsv_h2.gz -g 1.4e+09 -f BED --nomodel --extsize 200 --shift -100 -n H2AZ_KD --outdir ...../Peak_calling" ', 
                   wait = TRUE,  ignore.stderr = FALSE,  ignore.stdout = FALSE)

peak.path_h2 <- "..../H2AZ_KD_peaks.narrowPeak"
h2.peaks <- rtracklayer::import(peak.path_h2, format = "narrowPeak")
h2.peaks <- keepStandardChromosomes(h2.peaks, pruning.mode = 'coarse')

macs_count.h2 <- FeatureMatrix(fragments = Fragments(H2.data),
                               features = h2.peaks,
                               cells = colnames(H2.data))
H2.data[['peaks']] <- CreateChromatinAssay(
  counts = macs_count.h2,
  sep = c(":", "-"),
  genome = genome(genome.h2),
  fragments = fragpath_h2,
  min.cells = 10,
  annotation = genome.h2
)

#QC
DefaultAssay(H2.data) <- "peaks"
# compute nucleosome signal score per cell
H2.data <- NucleosomeSignal(object = H2.data)

# compute TSS enrichment score per cell
H2.data <- TSSEnrichment(object = H2.data, fast = FALSE)
DensityScatter(H2.data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

saveRDS(H2.data, file = "sox10ncc_h2az_kd.rds")

