library(Seurat)
library(scales)
library(ggplot2)
library(readr)
library(phytools)
library(treeio)
library(ggtree)
library(ape)         # For phylogenetic tree tools
library(phangorn)
#load these datasets
#you'll need the full seurat dataset
pbmc.combined = readRDS("~/SC Experiments/20230622_LN_GC/20230701_clean_demux/20230715_full_merged_enclone_gc.rds")
#and the one with the cleaned cells and plasma cells removed
pbmc.combined.clean = readRDS("~/SC Experiments/20230622_LN_GC/20230701_clean_demux/20230915_merged_demux_clean_plasma_rm_enclone_gc.rds")
#
trees = read.tree("~/SC Experiments/20230622_LN_GC/d1_wt_1.nwk.txt")
pbmc.sub = subset(pbmc.combined,subset = sample == 'D1_WT_1')
enclone <- read_csv("~/SC Experiments/20230622_LN_GC/enclone_complete_pcell/enclone_d1_wt_1_complete_pcell_melt.csv")
#specify sample of interest
sample = 'D1_WT_1'
#specify which zipcodes to plot 1,2,3,4,5,6
#Unecessary for GC plotting 
barcodes.to.include = c(1,2,3,4,5,6)
barcodes.to.include = paste('Barcode-',barcodes.to.include,sep = '')
#specify which clonotypes you want:
clonotype = 1

#the rest 
enclone.sub = subset(enclone,subset = group_id == clonotype)
#load tree
test.tree = trees[[clonotype]]
#initialize the data frame df which has all tip information
df = data.frame(label = unique(enclone.sub$exact_subclonotype_id))
#generate the barcode list
barcode.list <- vector("list", length(df$label))
names(barcode.list) <- df$label
for (name in df$label) {
  barcode.list[[name]] <- enclone.sub$barcodes[enclone.sub$exact_subclonotype_id == name]
}
n = c()
GC.call = rep('ambig',length(df$label))
ZC.call = rep('ambig',length(df$label))
for (i in 1:length(df$label)) {
  barcodes.temp = barcode.list[[i]]
  #ZC.temp = pbmc.sub$ZC.call[match(barcodes.temp,colnames(pbmc.sub))]
  GC.temp = pbmc.sub$GC.call[match(barcodes.temp,colnames(pbmc.sub))]
  GC.temp[!(barcodes.temp %in% colnames(pbmc.combined.clean))]<-'ambig'
  n[i] = sum(GC.temp != 'ambig')
  GC.temp <- GC.temp[GC.temp != "ambig"]
  # Step 2: get most frequent value
  if(sum(GC.temp != 'ambig')>0){
    GC.call[i] <- names(which.max(table(GC.temp)))
    
  }
  ZC.temp = pbmc.sub$ZC.call[match(barcodes.temp,colnames(pbmc.sub))]
  ZC.temp[!(ZC.temp %in% barcodes.to.include)]<-'ambig'
  ZC.temp[!(barcodes.temp %in% colnames(pbmc.combined.clean))]<-'ambig'
  ZC.temp <- ZC.temp[ZC.temp != "ambig"]
  # Step 2: get most frequent value
  if(sum(ZC.temp != 'ambig')>0){
    ZC.call[i] <- names(which.max(table(ZC.temp)))
    
  }
}
#add isotype information
df$isotype = enclone.sub$const1[match(df$label,enclone.sub$exact_subclonotype_id)]
df$label.isotype = paste(df$label,df$isotype,sep = '-')
df$n = n
df$GC.call = GC.call
df$ZC.call = ZC.call

#df$ZC.call = factor(df$ZC.call,levels = c('ambig','Barcode-1','Barcode-2','Barcode-3','Barcode-4','Barcode-5','Barcode-6'))
df$label = as.character(df$label)
tree = full_join(test.tree, df , by='label')
cols = c('gray',hue_pal()(6))
indices.to.drop = tree@phylo$tip.label[tree@extraInfo$GC.call == 'ambig']
tree.2 = drop.tip(tree,tip = indices.to.drop)

#get ancestral node states
GC.calls <- tree.2@extraInfo$GC.call[!is.na(tree.2@extraInfo$GC.call)] %>% as.factor()
GC.calls <- droplevels(GC.calls)     # remove unused levels
df <- as.data.frame((GC.calls))
df = as.data.frame(t(df))
colnames(df)<-tree.2@phylo$tip.label
row.names(df)<-1
dat <- phyDat(df, type = "USER", levels = levels(GC.calls))
tr_phylo <- as.phylo(tree.2) 
minSteps <- parsimony(tr_phylo, dat, method = "fitch")
anc <- ancestral.pars(tr_phylo, dat, type = "ACCTRAN")
anc_mat <- as.matrix(anc)
#node_states <- apply(anc_mat, 1, function(x) names(which(x == max(x))))
#edges <- tr_phylo$edge
levels_vec <- attr(dat, "levels")

anc_sampled <- lapply(anc, function(x) {
  probs <- as.numeric(x)
  probs <- probs / sum(probs)   # normalize just in case
  
  sample(levels_vec, size = 1, prob = probs)
})

anc_states <- lapply(anc, function(x) levels_vec[which(x == 1)])
n_tips <- length(tr_phylo$tip.label)
for(i in (n_tips+1):length(tree.2@extraInfo$GC.call)){
  tree.2@extraInfo$GC.call[i]<-anc_sampled[[i]]
}
#plot the tree
#all ancestral nodes will be colored based on inferred GC by parsimony
#y axis labels are the exact subclonotype designation
#tip/node labels are internal numbers for keeping track of each node, meaningless 

t = ggtree(tree.2,aes(color=GC.call),size=1) + geom_tiplab(linetype = 'dashed',as_ylab=TRUE,offset = .6, hjust = .5,linesize = 1) + geom_treescale()+
  geom_tippoint(aes(size = n,color = GC.call)) + geom_nodepoint(aes(color = GC.call),size=3,shape = 1,stroke=2)+
  geom_text(aes(label = node),color = 'black', hjust = -0.3)+ scale_size_area()+
  theme(legend.position = "right") +scale_color_manual(values =c('ambig'='gray',
                                                                 'GC1' = cols[2],
                                                                 'GC2' = cols[3],
                                                                 'GC3' = cols[4],
                                                                 'GC4' = cols[5],
                                                                 'GC5' = cols[6],
                                                                 'GC6' = cols[7]))+ggtitle(paste(paste(sample,'clonotype'),as.character(clonotype)))+
  theme_bw()


t
#plot the tree w out the node labels

t = ggtree(tree.2,aes(color=GC.call),size=1) + geom_tiplab(linetype = 'dashed',as_ylab=TRUE,offset = .6, hjust = .5,linesize = 1) + geom_treescale()+
  geom_tippoint(aes(size = n,color = GC.call)) + geom_nodepoint(aes(color = GC.call),size=3,shape = 1,stroke=2)+
  scale_size_area()+
  theme(legend.position = "right") +scale_color_manual(values =c('ambig'='gray',
                                                                 'GC1' = cols[2],
                                                                 'GC2' = cols[3],
                                                                 'GC3' = cols[4],
                                                                 'GC4' = cols[5],
                                                                 'GC5' = cols[6],
                                                                 'GC6' = cols[7]))+ggtitle(paste(paste(sample,'clonotype'),as.character(clonotype)))+
  theme_bw()


t

