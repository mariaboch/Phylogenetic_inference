#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ape)
library(ggplot2)

# paths

MAINDIR <- getwd()
figurespath <- paste(MAINDIR, "figures", sep="/")
consensuspath <- paste(MAINDIR, "data/consensus", sep="/")
supertreepath <- paste(MAINDIR, "data/supertree", sep="/")

# read trees

input_path <- paste(MAINDIR, "Mycobacterium_species_tree.nwk", sep="/")
species_tree <- read.tree(input_path)

input_path <- paste(supertreepath, "super_tree_0.nwk", sep="/")
basic_super <- read.tree(input_path)

input_path <- paste(consensuspath, "majority_consensus.nwk", sep="/")
basic_consensus <- read.tree(input_path)

if ("-b" %in% args) {
  input_path <- paste(supertreepath, "super_tree_bootstrap_0.nwk", sep="/")
  bootstrap_super <- read.tree(input_path)
  
  input_path <- paste(consensuspath, "majority_consensus_bootstrap.nwk", sep="/")
  bootstrap_consensus <- read.tree(input_path)
}

if ("-p" %in% args) {
  input_path <- paste(supertreepath, "super_tree_paralogs_0.nwk", sep="/")
  paralogs_super <- read.tree(input_path)
}

# plotting

# comparisons

#####################################################################
#####################################################################

out_path <- paste(figurespath, "Species_vs_basic_super.png", sep="/")

png(file=out_path,
    width = 20, height = 10, # Width and height in inches,
    unit = "in",
    res = 180,
    par(mar=c(0,0,0,0)),
    )

comparePhylo(species_tree, basic_super,
             plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")

dev.off()

#####################################################################

out_path <- paste(figurespath, "Species_vs_basic_consensus.png", sep="/")

png(file=out_path,
    width = 20, height = 10, # Width and height in inches,
    unit = "in",
    res = 180,
    par(mar=c(0,0,0,0)),
)

comparePhylo(species_tree, basic_consensus, force.rooted = TRUE,
             plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")

dev.off()


#####################################################################
# Super vs consensus

out_path <- paste(figurespath, "Basic_super_vs_basic_consensus.png", sep="/")

png(file=out_path,
    width = 20, height = 10, # Width and height in inches,
    unit = "in",
    res = 180,
    par(mar=c(0,0,0,0)),
)

comparePhylo(basic_super, basic_consensus, force.rooted = TRUE,
             plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")

dev.off()

#####################################################################
#####################################################################

if ("-b" %in% args) {
  
  #####################################################################

  out_path <- paste(figurespath, "Species_vs_bootstrap_super.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(species_tree, bootstrap_super,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################
  
  out_path <- paste(figurespath, "Species_vs_bootstrap_consensus.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(species_tree, bootstrap_consensus, force.rooted = TRUE,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################
  
  # Super vs consensus
  
  out_path <- paste(figurespath, "Bootstrap_super_vs_bootstrap_consensus.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(bootstrap_super, bootstrap_consensus, force.rooted = TRUE,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################
  
  # Consensus trees
  
  out_path <- paste(figurespath, "Basic_consensus_vs_bootstrap_consensus.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(basic_consensus, bootstrap_consensus,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  
  #####################################################################
  
  # Super
  
  out_path <- paste(figurespath, "Basic_super_vs_bootstrap_super.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(basic_super, bootstrap_super,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################

}

#####################################################################
#####################################################################

if ("-p" %in% args) {
  
  #####################################################################
  
  out_path <- paste(figurespath, "Species_vs_paralogs_super.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(species_tree, paralogs_super,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################
  # Super
  
  out_path <- paste(figurespath, "Basic_super_vs_paralogs_super.png", sep="/")
  
  png(file=out_path,
      width = 20, height = 10, # Width and height in inches,
      unit = "in",
      res = 180,
      par(mar=c(0,0,0,0)),
  )
  
  comparePhylo(basic_super, paralogs_super,
               plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
  
  dev.off()
  
  #####################################################################
  
  if ("-b" %in% args) {
  
    # Super
    
    out_path <- paste(figurespath, "Bootstrap_super_vs_paralogs_super.png", sep="/")
    
    png(file=out_path,
        width = 20, height = 10, # Width and height in inches,
        unit = "in",
        res = 180,
        par(mar=c(0,0,0,0)),
    )
    
    comparePhylo(bootstrap_super, paralogs_super,
                 plot = TRUE, label.offset = 0.8, align.tip.label=TRUE, lab4ut="axial", location="right")
    
    dev.off()
    
  #####################################################################
  
  }

}

#####################################################################
#####################################################################
