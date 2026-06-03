H2.combined <- readRDS("..../Integrated_H2.rds")
DefaultAssay(H2.combined) <- "SCT"
DimPlot(H2.combined, reduction = "umap.rna.h2", split.by = "orig.ident", label = T)
DotPlot(H2.combined, features = c("twist1a", "grem2b", "col11a1b", "lamc3", "foxd3", "crestin", "sox10", "zeb2a", "aox5", "paics", "dct", "mitfa", "tyr", "ltk", "alx4a", "gfap", "pou3f1", "her12", "elavl3", "elavl4", "oc90", "epcam", "neurod1", "musk", "tpma")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="viridis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() + 
  coord_flip() + ggtitle("Markers for H2A.Z Cluster Annotation")

#RNA Derived
Idents(H2.combined) <- H2.combined$integrated_snn_res.0.5
celltype <- rep(NA, length = ncol(H2.combined))

celltype[which(Idents(H2.combined) %in% c(0,10))] <- 'twist1a+ NCC'
celltype[which(Idents(H2.combined) %in% c(1))] <- 'mesenchymal'
celltype[which(Idents(H2.combined) %in% c(2))] <- 'unknown'
celltype[which(Idents(H2.combined) %in% c(4))] <- 'MX+'
celltype[which(Idents(H2.combined) %in% c(9))] <- 'pigment progenitor'
celltype[which(Idents(H2.combined) %in% c(8))] <- 'foxd3+ NCC'
celltype[which(Idents(H2.combined) %in% c(7))] <- 'melanoblast'
celltype[which(Idents(H2.combined) %in% c(13))] <- 'MI+'
celltype[which(Idents(H2.combined) %in% c(3))] <- 'neuronal'
celltype[which(Idents(H2.combined) %in% c(6,12))] <- 'otic'
celltype[which(Idents(H2.combined) %in% c(5))] <- 'glial/neural'
celltype[which(Idents(H2.combined) %in% c(11))] <- 'muscle progenitor'
celltype <- factor(celltype, 
                   levels = c('twist1a+ NCC', 'mesenchymal', 'foxd3+ NCC', 'pigment progenitor', 'MX+', 
                              'melanoblast', 'MI+', 'unknown', 'glial/neural',
                              'neuronal', 'otic', 'muscle progenitor'), 
                   ordered = T)

H2.combined$celltype.rna_n <- celltype
color_layers= c('mesenchymal' ='#efd129',
                'twist1a+ NCC'= '#e24041',
                'neuronal' ='#182953',
                'glial/neural' ='#277ea6',
                'unknown'= '#3ac4e7',
                'otic' ='#ab93c6',
                'foxd3+ NCC' = '#b06500',
                'pigment progenitor' = '#808000',
                'MX+'= '#7ac143',
                'melanoblast'= '#40733e',
                'MI+'='#ff9900',
                'muscle progenitor' = '#f64a8a'
)
Idents(H2.combined) <- H2.combined$celltype.rna_n
p_rna <- DimPlot(H2.combined, reduction = "umap.rna.h2", split.by = "orig.ident", pt.size = 1.2, cols = color_layers)
p_rna
saveRDS(H2.combined, file = "E:/Transit/Single_Cell+ATAC/New_Integrative_Analysis/H2AZ_sample/Data/07.08.2025/Integrated_H2_annotated.rds")

