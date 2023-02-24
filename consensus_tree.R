#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ape)

MAINDIR <- getwd()

consensuspath <- paste(MAINDIR, "data/consensus", sep="/")

input_path <- paste(consensuspath, "trees_for_consensus.nwk", sep="/")
output_path <- paste(consensuspath, "majority_consensus.nwk", sep="/")

#majority consensus without bootstrapping
all_nj_trees <- read.tree(input_path, keep.multi = TRUE)
majority_consensus <- consensus(all_nj_trees, p = 0.5, check.labels = TRUE, rooted = FALSE)
write.tree(majority_consensus, file = output_path, append = FALSE,
           digits = 10, tree.names = FALSE)

if (args[1]=="-b") {
  
  input_path_bootstrap <- paste(consensuspath, "trees_for_consensus_bootstrap.nwk", sep="/")
  output_path_bootstrap <- paste(consensuspath, "majority_consensus_bootstrap.nwk", sep="/")
  
  #majority consensus with bootstrapping
  all_nj_trees <- read.tree(input_path_bootstrap, keep.multi = TRUE)
  majority_consensus <- consensus(all_nj_trees, p = 0.5, check.labels = TRUE, rooted = FALSE)
  write.tree(majority_consensus, file = output_path_bootstrap, append = FALSE,
             digits = 10, tree.names = FALSE)

}
