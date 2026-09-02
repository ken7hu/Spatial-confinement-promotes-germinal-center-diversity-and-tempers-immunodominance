
# Used in Extended FIgure 4E ----------------------------------------------


library(Seurat)
library(ggplot2)
library(magrittr)
# data processing, unecessary for plotting --------------------------------
#code for generating the venn diagrams of DE genes
#read in the seurat object
pbmc.combined = readRDS("~/SC Experiments/20230622_LN_GC/20230701_clean_demux/20230713_full_merged_enclone.rds")
Idents(pbmc.combined)<-pbmc.combined$ln
markers_d1 = FindMarkers(pbmc.combined,ident.1 = 'D1_WT',ident.2 = 'D1_KO',min.pct = 0.1,logfc.threshold = 0.2)
markers_d2 = FindMarkers(pbmc.combined,ident.1 = 'D2_WT',ident.2 = 'D2_KO',min.pct = 0.1,logfc.threshold = 0.2)
markers_d3 = FindMarkers(pbmc.combined,ident.1 = 'D3_WT',ident.2 = 'D3_KO',min.pct = 0.1,logfc.threshold = 0.2)
#For venn diagramming of DE genes
d3_wt = subset(markers_d1, subset = p_val<0.05) %>% subset(subset = avg_log2FC >0) %>% row.names()
d3_ko = subset(markers_d1, subset = p_val<0.05) %>% subset(subset = avg_log2FC <0) %>% row.names()

d2_wt = subset(markers_d2, subset = p_val<0.05) %>% subset(subset = avg_log2FC >0) %>% row.names()
d2_ko = subset(markers_d2, subset = p_val<0.05) %>% subset(subset = avg_log2FC <0) %>% row.names()

d1_wt = subset(markers_d3, subset = p_val<0.05) %>% subset(subset = avg_log2FC >0) %>% row.names()
d1_ko = subset(markers_d3, subset = p_val<0.05) %>% subset(subset = avg_log2FC <0) %>% row.names()

library(ggVennDiagram)
wt.genes = list(d1_wt,d2_wt,d3_wt)
ko.genes = list(d1_ko,d2_ko,d3_ko)
names(wt.genes) = c('D1','D2','D3')

# plotting begins here ----------------------------------------------------
ggVennDiagram(wt.genes, label_alpha = 0,set_size = 5)+theme(legend.position = 'None')+scale_fill_distiller(palette = "Reds")
ggVennDiagram(ko.genes, label_alpha = 0,set_size = 5)+theme(legend.position = 'None')+scale_fill_distiller(palette = "Blues")

