# Phylogenetic_inference

## Final assignment for comparative genomics classes
## author: Maria Bochenek

This pipeline infers phylogenetic tree for collection of species genomes names given in a similar format as `filename`. It downloads genomes, retrieves genes, then clusters the genes into homology groups. Only "1-1" homology groups are selected for further inference. Then program computes multiple sequence alignment, infers NJ trees for each group following two cases: 
* without bootstraping and 
* with bootstraping alignments and selecting for further analysis only trees with mean bootstrap support $>=$ 75.

Finally program infers species tree using consensus and supertree methods for both sets of NJ trees (without bootstraping, with bootstraping).

## Required software and packages

* python: Biopython
* [MMSeq2](https://github.com/soedinglab/MMseqs2)
* [ClustalW](http://www.clustal.org/clustal2/)
* [fasturec](http://bioputer.mimuw.edu.pl/gorecki/fasturec/) and add its executable to your PATH or place executable in the same working directory as `phylogeny.py`
* R: ape

## Usage
`phylogeny.py`
