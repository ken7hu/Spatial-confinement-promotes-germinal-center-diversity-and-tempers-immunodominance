#read in the all_GC_overlaps_w_distances.rds file
total.df = readRDS(file = "~/SC Experiments/20230622_LN_GC/all_GC_overlaps_w_distances.rds")
library(RColorBrewer)
#assign reds to KO and blues to WT
reds = brewer.pal(n=9,name="Reds")
blues = brewer.pal(n=9,name="Blues")

#with samples colored differently
ggplot(data = total.df)+geom_point(aes(x=dist.value,y=overlap.value,color=sample),size=2)+ylim(0,3)+
  geom_hline(yintercept = 1,linetype='dashed',size=1)+geom_smooth(aes(x=dist.value,y=overlap.value,group=genotype,color = genotype),method='lm') +scale_color_manual(values =
                                                                                                                                                                       c(reds[8],reds[8],blues[9],blues[8],reds[7],blues[7],reds[6],reds[5],blues[6],blues[5],'red','blue'))+
  ylab('Overlap Enrichment Between GCs')+xlab('Distance between GCs (microns)')+theme(axis.title = element_text(size=14))+
  theme_bw()+theme(axis.title = element_text(size=16),axis.text = element_text(size=12))


#pooled by genotype
ggplot(data = total.df)+geom_point(aes(x=dist.value,y=overlap.value,color=genotype),size=2)+ylim(0,3)+
  geom_hline(yintercept = 1,linetype='dashed',size=1)+geom_smooth(aes(x=dist.value,y=overlap.value,group=genotype,color = genotype),method='lm') +scale_color_manual(values =
                                                                                                                                                                       c(reds[8],blues[8]))+
  ylab('Overlap Enrichment Between GCs')+xlab('Distance between GCs (microns)')+theme(axis.title = element_text(size=14))+
  theme_bw()+theme(axis.title = element_text(size=16),axis.text = element_text(size=12))
