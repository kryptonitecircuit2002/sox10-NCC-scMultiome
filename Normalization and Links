#RNA Assay normalization
data<- readRDS(".../preprocessed_data.rds")
data <-  SCTransform(data, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

#Peaks normalisation
DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = "q5")
data <- RunSVD(data)
seqlevelsStyle(data) <- "NCBI"
seqlevelsStype(BSgenome.Drerio.UCSC.danRer11) <- "NCBI"  #To ensure seqlevels match between peaks of object and supplied genome
data <- RegionStats(data, genome = BSgenome.Drerio.UCSC.danRer11)
data <- LinkPeaks(data,
                  peak.assay = "peaks",
                  expression.assay = "SCT",
                  pvalue_cutoff = 1,
                  score_cutoff = 0,
                  method = 'spearman')

## Include adjusted p-values information into links
links <- data[['peaks']]@links
pvalue_adjusted = p.adjust(links$pvalue, method = 'fdr')
links$pvalue.fdr = pvalue_adjusted
data[['peaks']]@links <- links

##Filter links
links <- links[intersect(which(links$pvalue.fdr < 0.1), which(abs(links$score) > 0.1)),]

# 2. seperating promoter and enhancer linked peaks
gene.ranges <- CollapseToLongestTranscript(ranges = Annotation(data)) ##CollapseToLongestTranscript is a custom defined function (see below)
peaks <- StringToGRanges(links$peak, sep = c("-", "-"))
gene.promoters <- promoters(gene.ranges, upstream = 500, downstream = 500)
hits <- findOverlaps(query = peaks, subject = gene.promoters)
promoter.peaks <- unique(queryHits(hits))
links <- links[-promoter.peaks] 

save(links, file = "..../enhancer_links_filtered_allcells.RData")
save(gene.promoters, file = "..../promoterpeaks_allcells.RData")

##ATAC assay
DefaultAssay(data) <- "ATAC"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = "q5")
data <- RunSVD(data)
data <- RegionStats(data, genome = BSgenome.Drerio.UCSC.danRer11)
data <- LinkPeaks(data,
                  peak.assay = "ATAC",
                  expression.assay = "SCT",
                  pvalue_cutoff = 1,
                  score_cutoff = 0,
                  method = 'spearman')
## Include adjusted p-values information into links
ATAC_links <- data[['ATAC']]@links
pvalue_adjusted = p.adjust(ATAC_links$pvalue, method = 'fdr')
ATAC_links$pvalue.fdr = pvalue_adjusted
data[['ATAC']]@links <- ATAC_links

##Filter links
ATAC_links <- ATAC_links[intersect(which(ATAC_links$pvalue.fdr < 0.1), which(abs(ATAC_links$score) > 0.1)),]

# 2. seperating promoter and enhancer linked peaks
ATAC_peaks <- StringToGRanges(ATAC_links$peak, sep = c("-", "-"))
ATAC_hits <- findOverlaps(query = ATAC_peaks, subject = gene.promoters)
promoter.peaks_ATACassay <- unique(queryHits(ATAC_hits))
ATAC_links <- ATAC_links[-promoter.peaks_ATACassay]

save(ATAC_links, file = "..../peaksassay_enhancer_links_filtered_allcells.RData")
save(gene.promoters, file = "..../promoterpeaks_allcells.RData")

# PCA, clustering, and UMAP visualization
DefaultAssay(data) <- "SCT"
data <- RunPCA(data, verbose = FALSE, reduction.name = "pca")
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.7)
data <- RunUMAP(data, reduction = "pca", reduction.name = "umap.rna", dims = 1:30, verbose = FALSE, spread = 0.32, min.dist = 0.35)
DimPlot(data, reduction = "umap.rna", label = TRUE)

#Save RData
saveRDS(data, file = ".../data_linked.rds")


CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ]
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
