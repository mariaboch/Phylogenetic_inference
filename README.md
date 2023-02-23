# Phylogenetic_inference

## Final assignment for comparative genomics classes
## author: Maria Bochenek

This pipeline infers phylogenetic tree for collection of species genomes names given in a similar format as `Mycobacterium_list.txt`. It downloads genomes, retrieves genes, then clusters the genes into homology groups. Only "$1-1$" homology groups are selected for further inference. Then program computes multiple sequence alignment, infers NJ trees for each group following two cases: 
* without bootstraping and 
* with bootstraping alignments and selecting for further analysis only trees with mean bootstrap support $\geq$ 75.

Finally program infers species tree using consensus and supertree methods for both sets of NJ trees (without bootstraping, with bootstraping).

## Required software and packages

* python: Biopython
* R: ape
* [MMSeq2](https://github.com/soedinglab/MMseqs2)
* [ClustalW](http://www.clustal.org/clustal2/)
* [fasturec](http://bioputer.mimuw.edu.pl/gorecki/fasturec/) and add its executable to your PATH or place executable in the same working directory as `phylogeny.py`

## Usage
`python3 phylogeny.py [-h] -n genomelistfile [--paralogs] [--bootstap] [-thr] [-e]`

where
* `-h, --help` show help message and exit
* `-n genomelistfile path` to a file containing species' genome names relative to main directory
* `--paralogs` calculate supertree including clusters with paralogs
* `--bootstap` calculate bootstrap supports for gene trees.
* `-bn` number of bootstrap times default=50.
* `-thr, --mean_support_thr` Mean bootstrap support threshold. Trees with mean bootstrap support below threshold will be discarded. Default value is 75 (75%).
* `-e, --email` E-mail address (for Entrez download)
