# Usage: run command line Rscript midpoint_tree.R gene
# For example: Rscript scripts/midpoint_tree.R YOL164W

# Midpoint root a phylogeny with R:
library(ape)
library(phangorn)

args = commandArgs(T)
print(args[1])

treefile <- paste("./intermediate/",args[1],".treefile",sep="")
mytree <- read.tree(treefile)
midponit_mytree <- ladderize (midpoint(mytree))
write.tree (midponit_mytree, paste("./intermediate/",args[1],"_midpoint.tree",sep=""))
