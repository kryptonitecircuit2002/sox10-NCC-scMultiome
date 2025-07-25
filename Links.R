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
library(data.table)

set.seed(1234)
#Read Data
data<- readRDS(".../data_clustered.rds")

#Peaks Assay 
DefaultAssay(data) <- "peaks"
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
