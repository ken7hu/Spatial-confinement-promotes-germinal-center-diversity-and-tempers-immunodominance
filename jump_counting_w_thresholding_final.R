library(ggpubr)
library(readxl)
library(readr)
# For plotting the number of GC jumps when filtering out jumps that occur < 1 evo distance from root node -----------------------------------------
#load in the compiled dataframe depending on cutoff of evo distance you want (0.5 or 1) from compiled_0.5_cutoff.csv

#or compiled_1_cutoff.csv
compiled = read.csv("~/SC Experiments/20230622_LN_GC/jump_records/compiled_1_cutoff.csv")


# actual plotting code ----------------------------------------------------
#Not in the ppt
#plot for each sample plot all clonotypes integrated jump distances between GCs (microns) normalized to nNodes but w jumps that are too close to root filtered out
ggboxplot(
  compiled, x = "sample", y = "dist_integrated_norm_nnode",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6, outlier.shape = NA
) +
  geom_jitter(
    aes(x = sample, y = dist_integrated_norm_nnode, color = genotype),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15),
    size = 2, alpha = 0.8
  ) +
  labs(x = "Sample", y = "jump distances/nNodes",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# Slide 8
# plot for each sample the number of jumps between GC normalized to nNodes but w jumps that are too close to root filtered out (Slide 8)
ggboxplot(
  compiled, x = "sample", y = "norm.jumps.cutoff",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6, outlier.shape = NA
) +
  geom_jitter(
    aes(x = sample, y = norm.jumps.cutoff, color = genotype),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15),
    size = 2, alpha = 0.8
  ) +
  labs(x = "Sample", y = "jumps/nNodes",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#not in the ppt
# plot all integrated jump distances between GC w jumps that are too close to root filtered out but pooled by genotype ------------------------
compiled$genotype <- factor(compiled$genotype)
my_comparisons <- combn(levels(compiled$genotype), 2, simplify = FALSE)
ggboxplot(
  compiled, x = "genotype", y = "dist_integrated_norm_nnode",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6, outlier.shape = NA,
  add = "jitter",
  add.params = list(size = 2, alpha = 0.8,
                    position = position_jitter(width = 0.15))
) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",      # or "p.signif"
    hide.ns = TRUE
  ) +
  labs(x = "Genotype", y = "Jump Dist./nNode",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  rotate_x_text(45) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) 
# Slide 7
# plot pooled avg jump number normalzied to node (Slide 7)------------------------
compiled$genotype <- factor(compiled$genotype)

my_comparisons <- combn(levels(compiled$genotype), 2, simplify = FALSE)
ggboxplot(
  compiled, x = "genotype", y = "norm.jumps.cutoff",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6, outlier.shape = NA,
  add = "jitter",
  add.params = list(size = 2, alpha = 0.8,
                    position = position_jitter(width = 0.15))
) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",      # or "p.signif"
    hide.ns = TRUE
  ) +
  labs(x = "Genotype", y = "Jumps/nNode",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  rotate_x_text(45) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) 



# Prep for plotting on Slides 9-11 for plotting pooled evo distances from root for every jump pooled by genotype-------------------------
#set the cutoff (if you dont want any filtering, set cutoff to -1)
cutoff = -1
#read in the file of jump distances called jump.evo.dist.csv
jump.evo.dist = read.csv("~/SC Experiments/20230622_LN_GC/jump.evo.dist.csv")

#Slide 9
#for plotting evo distance for each GC jump but split out by samples
ggboxplot(
  subset(jump.evo.dist,subset = dist_from_root>1), x = "sample", y = "dist_from_root",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6, outlier.shape = NA
) +
  geom_jitter(
    aes(x = sample, y = dist_from_root, color = genotype),
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15),
    size = 2, alpha = 0.8
  ) +
  labs(x = "Sample", y = "Evo dist from root",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Slides 10 (cutoff set to -1)and 11 (cutoff set to 1) -------------------------------------------------------------
#for plotting evo distance for each GC jump but split out by genotypes pooling across samples

jump.evo.dist$genotype <- factor(jump.evo.dist$genotype)

my_comparisons <- combn(levels(jump.evo.dist$genotype), 2, simplify = FALSE)
ggboxplot(
  subset(jump.evo.dist,subset = dist_from_root>cutoff), x = "genotype", y = "dist_from_root",
  color = "genotype", fill = "genotype",
  alpha = 0.5, width = 0.6,
  add = "jitter",
  add.params = list(size = 2, alpha = 0.8,
                    position = position_jitter(width = 0.15))
) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",      # or "p.signif"
    hide.ns = TRUE
  ) +
  labs(x = "Genotype", y = "Evo dist from root",
       color = "Genotype", fill = "Genotype") +
  theme_pubr() +
  rotate_x_text(45) 
