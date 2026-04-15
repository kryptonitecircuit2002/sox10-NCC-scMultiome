lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

data <- readRDS(".../data_annotated.rds")
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
data <- AddMotifs(
  object = data,
  genome = BSgenome.Drerio.UCSC.danRer11,
  pfm = pfm
)
#Find differentially accessible regions
DefaultAssay(data) <- "peaks"
dar <- FindAllMarkers(
  object = data,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
top.dar <- rownames(dar[dar$p_val < 0.005 & dar$pct.1 > 0.2, ]) #filter dars

#Perform TF Motif Enrichment
enriched.motifs <- FindMotifs(
  object = data,
  features = top.dar)
enriched.motifs <- enriched.motifs[enriched.motifs$p.adjust < 0.05, ]

MotifPlot(
  object = data,
  motifs = head(rownames(enriched.motifs))
)

saveRDS(".../data_motifenriched.rds")
