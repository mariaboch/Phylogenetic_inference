#!/usr/bin/env python3

import os
import subprocess
import argparse
import matplotlib.pyplot as plt

# Bio
import Bio
from Bio import Phylo

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Builds consensus tree and supertree for a collection of species' genome names.")
    parser.add_argument("-n", required=True,
                        help="path to a file containing species' genome names relative to main working directory")
    parser.add_argument('-p', '--paralogs', dest="p", action='store_true', default=False,
                        help="Visualize paralogs tree")
    parser.add_argument('-b', '--bootstrap', dest="b", action='store_true', default=False,
                        help="Visualize bootstrap trees")
    args = parser.parse_args()
    return args

args = parse_arguments()

MAINDIR = os.getcwd()
DATADIR = MAINDIR+"/data"

GENOMELISTFILE = MAINDIR + f"/{args.n}"

datapath = DATADIR
supertreepath = datapath + "/supertree"
consensustreepath = datapath + "/consensus"

figurespath = MAINDIR+"/figures"
if not os.path.exists(figurespath):
    os.makedirs(figurespath)

def map_genome_id_to_name(file):
    mapping = {}
    with open(file, "r") as f:
        lines = [l.strip().split("_") for l in f.readlines()]
    f.close()
    for genome in lines:
        genome_id = "_".join(genome[-2:])
        if genome[-5]=='complete':
            genome_name = "M." + " ".join(genome[1:-5])
        else:
            genome_name = "M." + " ".join(genome[1:-2])
        mapping[genome_id] = genome_name
    return mapping

genome_id_to_name = map_genome_id_to_name(GENOMELISTFILE)

##########################################################################

# Species tree
tree = Phylo.read(f"{MAINDIR}/Mycobacterium_species_tree.nwk", "newick")
tree.ladderize()
for t in tree.get_terminals():
    t.name = genome_id_to_name[t.name]
fig = plt.figure(figsize=(10, 12), dpi=150)
axes = fig.add_subplot(1, 1, 1)
plt.title("Mycobacterium species tree", fontsize=18)
plt.xlabel("branch length", fontsize=15)
plt.ylabel("taxa", fontsize=15)
Phylo.draw(tree, axes=axes, do_show=False)
plt.savefig(f'{figurespath}/Mycobacterium_species_tree.png')

##########################################################################

#Supertree 
tree = Phylo.read(f"{supertreepath}/super_tree_0.nwk", "newick")
tree.ladderize()
for t in tree.get_terminals():
    t.name = genome_id_to_name[t.name]
fig = plt.figure(figsize=(10, 12), dpi=150)
axes = fig.add_subplot(1, 1, 1)
plt.title("Supertree", fontsize=18)
plt.xlabel("branch length", fontsize=15)
plt.ylabel("taxa", fontsize=15)
Phylo.draw(tree, axes=axes, do_show=False)
plt.savefig(f'{figurespath}/Basic_supertree.png')

##########################################################################

# Consensus
tree = Phylo.read(f"{consensustreepath}/majority_consensus.nwk", "newick")

tree.ladderize()
for t in tree.get_nonterminals():
    t.confidence = round(t.confidence, 3)
for t in tree.get_terminals():
    t.name = genome_id_to_name[t.name]
fig = plt.figure(figsize=(10, 12), dpi=150)
axes = fig.add_subplot(1, 1, 1)
plt.title("Consensus tree", fontsize=18)
plt.xlabel("branch length", fontsize=15)
plt.ylabel("taxa", fontsize=15)
Phylo.draw(tree, axes=axes, do_show=False)
plt.savefig(f'{figurespath}/Basic_consensus.png')

##########################################################################

if args.b:
    
    #Supertree
    tree = Phylo.read(f"{supertreepath}/super_tree_bootstrap_0.nwk", "newick")
    tree.ladderize()
    for t in tree.get_terminals():
        t.name = genome_id_to_name[t.name]
    fig = plt.figure(figsize=(10, 12), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    plt.title("Bootstrap supertree", fontsize=18)
    plt.xlabel("branch length", fontsize=15)
    plt.ylabel("taxa", fontsize=15)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(f'{figurespath}/Bootstrap_supertree.png')

    #Consensus
    tree = Phylo.read(f"{consensustreepath}/majority_consensus_bootstrap.nwk", "newick")
    tree.ladderize()
    for t in tree.get_nonterminals():
        t.confidence = round(t.confidence, 3)
    Phylo.write(tree, f"{consensustreepath}/majority_consensus_bootstrap_names.nwk", "newick")
    for t in tree.get_terminals():
        t.name = genome_id_to_name[t.name]
    fig = plt.figure(figsize=(10, 12), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    plt.title("Bootstrap consensus tree", fontsize=18)
    plt.xlabel("branch length", fontsize=15)
    plt.ylabel("taxa", fontsize=15)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(f'{figurespath}/Bootstrap_consensus.png')
##########################################################################

if args.p:

    tree = Phylo.read(f"{supertreepath}/super_tree_paralogs_0.nwk", "newick")
    tree.ladderize()
    for t in tree.get_terminals():
        t.name = genome_id_to_name[t.name]
    fig = plt.figure(figsize=(10, 12), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    plt.title("Paralogs supertree", fontsize=18)
    plt.xlabel("branch length", fontsize=15)
    plt.ylabel("taxa", fontsize=15)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(f'{figurespath}/Paralogs_supertree.png')

##########################################################################

if args.b and args.p:
    subprocess.call('/usr/bin/Rscript compare_trees.R -b -p', shell=True)
else:
    if args.b:
        subprocess.call('/usr/bin/Rscript compare_trees.R -b', shell=True)
    elif args.p:
        subprocess.call('/usr/bin/Rscript compare_trees.R -p', shell=True)
    else:
        subprocess.call('/usr/bin/Rscript compare_trees.R', shell=True)
