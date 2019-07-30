# check if libraries are installed 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggtree"))
  BiocManager::install("ggtree")

# load libraries 
require("ggtree")
library("ggplot2")
require("cowplot")
library("RColorBrewer")

# Inputs 
cwd <- "/Users/alvin/Dropbox/github_repo/Phydelity/manuscript/FLU" # working directory 
newick_tree_file <- "raxml/phydelity/zero-branch-length-collapsed_RAxML_bestTree.nwk" # path to newick tree file relative to working directory 
meta_df_file <- "phydelity_result_meta.csv"
outfname <- "supplementary_figure_5.pdf"

setwd(cwd)
dd = read.csv(meta_df_file, sep = ",", 
              stringsAsFactors = T, 
              na.strings = "")
# change numeric labels as characters 
dd$Phydelity = as.character(dd$Phydelity)
dd$tp = as.character(dd$tp)
dd$Household = as.character(dd$Household)

# take all R colors from graphical devices (removing shades of greys - 433 colours in total)
# and sample from them
set.seed(668)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols = sample(color, length(unique(na.omit(dd$tp)))+length(unique(na.omit(dd$Phydelity)))+length(unique(na.omit(dd$Household))))

# parse newick tree file 
tree <- read.newick(newick_tree_file) 
p <- ggtree(tree, aes(x, y)) + 
  geom_tree() + 
  theme_tree() + 
  geom_treescale(x=0, y=-2, fontsize = 2.8, width=0.001) 

# highlight tip labels of true transmission pairs
p <- p %<+% dd + geom_tiplab(aes(fill=tp),
                             color = "black", # color for label font
                             geom = "label",  # labels not text
                             label.padding = unit(0.1, "lines"), # amount of padding around the labels
                             label.size = 0, 
                             size=2.3) 

rownames(dd) <- dd$leaf # make leaf row as index --> required for gheatmap
p <- gheatmap(p, dd[c("Phydelity", "Household")], offset=0.001, width=.2, 
              font.size=2.8, colnames_angle=0, hjust=0.5, 
              colnames_offset_y=-2) + theme(legend.position = "none") + scale_fill_manual(values = cols)
ggsave(outfname, device = "pdf", width = 25, height = 50, units ="cm")
