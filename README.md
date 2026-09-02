# Spatial-confinement-promotes-germinal-center-diversity-and-tempers-immunodominance
Code used to analyze scRNA-Seq data in "Spatial confinement promotes germinal center diversity and tempers immunodominance"


tree_ancestornodes_colored_evodist.R: used for plotting phylogenetic trees with tips colored by GC assignment and ancestor nodes colored by inferred GC assignment according to ACCTRAN as seen in Fig. 
DE_venndiagram.R: used for Extended Figure 4E to plot overlapping DE genes
jump_counting_w_thresholding_final.R: used for the boxplot comparisons of number of GC jumps in clonotypes in WT and KO samples but w removal of jumps that occur too close in evolutionary distance to the ancestor node. Used in Extended Figure 4G
gc_sharing_distances_final.R: used for plotting overlap index by distance between GCs in Figure 4D
