# Usage: run command line Rscript midpoint_tree.R gene
# For example: Rscript midpoint_tree.R YFR055W

# Midpoint root a phylogeny with R:
library(ape)
library(phangorn)

# Remember to set your own directory here
setwd("~/Documents/Le/HGTdetect/case/saccharomyces_cerevisiae_S288c/tree_sce_S288c")

args = commandArgs(T)

print(args[1])

treefile <- paste("./fasta/",args[1],"/",args[1],".treefile",sep="")
mytree <- read.tree(treefile)
midponit_mytree <- ladderize (midpoint(mytree))
write.tree (midponit_mytree, paste("./fasta/",args[1],"/",args[1],"_midpoint.tree",sep=""))

# # And then generate the font and color for iTOL
# perl create_iTOL_config.pl ../fasta/YFR055W/YFR055W_midpoint.tree