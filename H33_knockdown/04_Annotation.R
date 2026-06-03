lapply(required_packages, library, character.only = TRUE)
set.seed(1234)

H3.combined <- readRDS("..../h3_integrated_linked.rds")
DimPlot(H3.combined, reduction = "umap.rna.h3", split.by = "orig.ident", label = TRUE, label.size = 4, repel = T)

#RNA Derived
Idents(H3.combined) <- H3.combined$SCT_snn_res.0.5
celltype <- rep(NA, length = ncol(H3.combined))

celltype[which(Idents(H3.combined) %in% c(1,11))] <- 'twist1a+ NCC'
celltype[which(Idents(H3.combined) %in% c(0))] <- 'mesenchymal'
celltype[which(Idents(H3.combined) %in% c(6))] <- 'foxd3+ NCC'
celltype[which(Idents(H3.combined) %in% c(5))] <- 'pigment progenitor'
celltype[which(Idents(H3.combined) %in% c(3))] <- 'MX+'
celltype[which(Idents(H3.combined) %in% c(10))] <- 'melanoblast'
celltype[which(Idents(H3.combined) %in% c(15))] <- 'MI+'
celltype[which(Idents(H3.combined) %in% c(4))] <- 'unknown'
celltype[which(Idents(H3.combined) %in% c(2))] <- 'glial/neural'
celltype[which(Idents(H3.combined) %in% c(7))] <- 'neuronal'
celltype[which(Idents(H3.combined) %in% c(8,13,16))] <- 'otic'
celltype[which(Idents(H3.combined) %in% c(14))] <- 'immature neurons'
celltype[which(Idents(H3.combined) %in% c(9,12))] <- 'muscle progenitor'
celltype <- factor(celltype, 
                          levels = c('twist1a+ NCC', 'mesenchymal', 'foxd3+ NCC', 'pigment progenitor', 'MX+', 
                                     'melanoblast', 'MI+', 'unknown', 'glial/neural',
                                     'neuronal', 'otic', 'immature neurons', 'muscle progenitor'), 
                          ordered = T)

H3.combined$celltype.rna_n <- celltype
Idents(H3.combined) <- H3.combined$celltype.rna_n
color_layers= c('mesenchymal' ='#efd129',
                'twist1a+ NCC'= '#e24041',
                'neuronal' ='#182953',
                'glial/neural' ='#277ea6',
                'immature neurons'= '#b22987',
                'unknown'= '#3ac4e7',
                'otic' ='#ab93c6',
                'foxd3+ NCC' = '#b06500',
                'pigment progenitor' = '#808000',
                'MX+'= '#7ac143',
                'melanoblast'= '#40733e',
                'MI+'='#ff9900',
                'muscle progenitor' = '#f64a8a'
)
p1 <- DimPlot(H3.combined, reduction = "umap.rna.h3", split.by = "orig.ident", cols = color_layers, pt.size = 1.1)
p1

#DotPlot
DefaultAssay(H3.combined.sct) <- "SCT"
DotPlot(H3.combined.sct, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "crestin", "sox10", "zeb2a", "aox5", "paics", "dct", "mitfa", "tyr", "ltk", "alx4a", "gfap", "pou3f1", "her12", "elavl3", "elavl4", "oc90", "epcam", "neurod1", "musk", "tpma")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() + 
  coord_flip() + ggtitle("Markers for Cluster Annotation")
saveRDS(H3.combined.sct, file = "..../H3_combined_annotated.rds")
