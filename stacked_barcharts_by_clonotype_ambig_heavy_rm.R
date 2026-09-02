#filter on clonotypes that have either IGL or IGK and IGH
library(readr)
library(ggplot2)
library(dplyr)
pbmc.combined.clean = readRDS("~/SC Experiments/20230622_LN_GC/20230701_clean_demux/20230915_merged_demux_clean_plasma_rm_enclone_gc.rds")
#specify which sample you want
sample.oi = 'D3_WT_2'
#the number of clonotypes to display
n.display = 30
#the enlone file should also match, change to match the sample.oi
enclone <- read_csv("~/SC Experiments/20230622_LN_GC/enclone_complete_pcell/enclone_d3_wt_2_complete_pcell_melt.csv")
#the rest should run by itself
pbmc.sub = subset(pbmc.combined.clean,subset = sample == sample.oi)
#remove all cells that ? in their enclone table column const1
keep = enclone$const1[match(colnames(pbmc.sub),enclone$barcodes)] != '?'
pbmc.sub = subset(pbmc.sub,cells = colnames(pbmc.sub)[keep])
####
pbmc.sub$GC.call = factor(pbmc.sub$GC.call,levels = c('GC1','GC2','GC3','GC4','GC5','GC6'))
clonotype.table =table(pbmc.sub$enclone.clone)[order(table(pbmc.sub$enclone.clone),decreasing = TRUE)]
clonotype.table = subset(clonotype.table,names(clonotype.table) != "None")
#remove any clonotypes with less than 2 chains
#assign k and l 
enclone$light.chain = rep("ambig",nrow(enclone))
enclone$light.chain[which(enclone$const1 == 'IGKC')]<-'kappa'
enclone$light.chain[which(enclone$const2 == 'IGKC')]<-'kappa'
enclone$light.chain[which(enclone$const3 == 'IGKC')]<-'kappa'

enclone$light.chain[which(enclone$const1 == 'IGLC1')]<-'lambda1'
enclone$light.chain[which(enclone$const2 == 'IGLC1')]<-'lambda1'
enclone$light.chain[which(enclone$const3 == 'IGLC1')]<-'lambda1'

enclone$light.chain[which(enclone$const1 == 'IGLC2')]<-'lambda2'
enclone$light.chain[which(enclone$const2 == 'IGLC2')]<-'lambda2'
enclone$light.chain[which(enclone$const3 == 'IGLC2')]<-'lambda2'

enclone$light.chain[which(enclone$const1 == 'IGLC3')]<-'lambda3'
enclone$light.chain[which(enclone$const2 == 'IGLC3')]<-'lambda3'
enclone$light.chain[which(enclone$const3 == 'IGLC3')]<-'lambda3'
#
good.clonotypes = enclone$group_id[enclone$nchains>1]
clonotype.table = subset(clonotype.table,names(clonotype.table) %in% good.clonotypes)
full.df = data.frame(Var1 = c(),Freq = c(),clonotype = c())
for(i in 1:n.display){
  clone = names(clonotype.table[i])
  #clone = names(clonotype.table[i+1])
  temp = table(pbmc.sub$GC.call[pbmc.sub$enclone.clone == clone]) %>% as.data.frame()
  temp$clonotype = rep(clone,nrow(temp))
  full.df = rbind(full.df,temp)
}
full.df$clonotype = factor(full.df$clonotype,levels = names(clonotype.table))
full.df$light.chain = enclone$light.chain[match(full.df$clonotype,enclone$group_id)]
full.df$heavy.chain = enclone$v_name1[match(full.df$clonotype,enclone$group_id)]
full.df$comb.chain = paste(full.df$heavy.chain,full.df$light.chain,sep = '/')
#get heights
df_bar_top <- full.df %>%
  group_by(clonotype) %>%
  summarise(
    total_y = sum(Freq)+12,
    comb.chain = comb.chain  # assuming 1 igv per sample
  )
#with the chain info annotated on top
ggplot(full.df, aes(fill=Var1, y=Freq, x=clonotype)) + 
  geom_bar(stat="identity")+theme_bw()+scale_y_continuous(expand = c(0.2, 0)) +
  geom_text(data = df_bar_top,aes(x = clonotype,y=total_y,label = comb.chain),inherit.aes = FALSE ,angle = 90)+
  theme(axis.title = element_blank(),axis.text = element_text(size=10),axis.text.x = element_text(angle = ,hjust=1,size=16),plot.margin = unit(c(1,1,1,1), "cm"))+
  ggtitle(sample.oi)

#without the chain info annotated
ggplot(full.df, aes(fill=Var1, y=Freq, x=clonotype)) + 
  geom_bar(stat="identity")+theme_bw()+scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title = element_blank(),axis.text = element_text(size=10),axis.text.x = element_text(angle = ,hjust=1,size=16),plot.margin = unit(c(1,1,1,1), "cm"))+
  ggtitle(sample.oi)


