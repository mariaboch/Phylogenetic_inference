#!/usr/bin/env python3

import os
import time
import re
from glob import glob
import subprocess
import argparse

# Bio
import Bio
from Bio import AlignIO, Entrez, GenBank
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.Consensus import bootstrap
from Bio.Phylo.Consensus import get_support
from Bio.Align.Applications import ClustalwCommandline

from ete3 import Tree

##########################################################################################################################
##########################################################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Builds consensus tree and supertree for a collection of species' genome names.")
    parser.add_argument("-n", required=True,
                        help="path to a file containing species' genome names relative to main working directory")
    parser.add_argument('-p', '--paralogs', dest="p", action='store_true', default=False,
                        help="Infere supertree using paralogs.")
    parser.add_argument('-b', '--bootstrap', dest="b", action='store_true', default=False,
                        help="Calculate bootstrap supports for gene trees.")
    parser.add_argument("-bn", default=50, type=float, required=False,
                        help="Number of bootstrap times default=50.")
    parser.add_argument("-thr", default=75, type=float, required=False,
                        help='Mean bootstrap support threshold (use with --bootstrap). Trees with mean bootstrap support below threshold will be discarded. Default value is 75 (75%).')
    parser.add_argument("-e", "--email", required=False, help="E-mail address (for Entrez download)",
                        default="m.bochenek@student.uw.edu.pl")
    args = parser.parse_args()
    return args

##########################################################################################################################
##########################################################################################################################

#Downloading the data 

def processGenomeList(file):
    genome_list = []
    newnames= open(file)
    for line in newnames.readlines():
        genomeId = line.split('_')
        genomeId = genomeId[-2] + '_' + genomeId[-1].strip('\n')
        genome_list.append(genomeId)
    return genome_list


def retrieve_genomes(EMAIL, GENOMELISTFILE):
    # Retrieve NCBI Complete Genome Records using Entez
    print("\nRetrieving NCBI complete genome records...")
    start = time.time()
    Entrez.email = EMAIL 

    genomeidlist = processGenomeList(GENOMELISTFILE)

    print("Number of genomes: ",str(len(genomeidlist)))
    for genomeid in genomeidlist:
            filename = "DB_"+genomeid+".gbk"
            handle = Entrez.efetch(db='nuccore', id=genomeid, rettype="gbwithparts", retmode="text")
            out_handle = open(filename, "w")
            for line in handle:
                    out_handle.write(line)
            out_handle.close()
            # print("Saved: ",filename)
    print("NCBI complete genome records retrieved in "+str(time.time()-start)+" seconds")


# Retrieving genes
def download_genes(genome_path, genes_path):
    cds_file = open(genes_path, "w")
    for file in os.listdir(genome_path):
        parser = GenBank.RecordParser()
        rc_file = open(genome_path + "/" + file)
        record = parser.parse(rc_file)
        rc_file.close()
        for feature in record.features:
            if feature.key == 'CDS':
                pid = ""
                trans = ""
                for qualifier in feature.qualifiers:
                    if qualifier.key == '/locus_tag=':
                        pid = qualifier.value[1:-1]
                        pid = pid.replace('"', '')
                    if qualifier.key == '/translation=':
                        trans = qualifier.value[1:-1]
                cds_file.write(">%s|%s" % (record.locus, pid + "\n"))
                cds_file.write("%s" % trans + "\n")
    cds_file.close()

##########################################################################################################################
##########################################################################################################################
# Reading clusters

def search_unwanted_aminoacids(string):
    # searching for aminoacids that are not in blossum62 matrix
    unwanted = "[J|U|O]"
    return re.search(unwanted, string)

def read_and_filter_clusters(genomeidlist, clusterpath, outclusterpath, filteredclusterpath, paralogsclusterpath=None):
    #creating directories
    if not os.path.exists(outclusterpath):
        os.makedirs(outclusterpath)
    if not os.path.exists(filteredclusterpath):
        os.makedirs(filteredclusterpath)

    #reading mmseq clusters
    with open(f"{clusterpath}/DB_clu_seq.fasta", "r") as f:
        lines = [l.strip() for l in f.readlines()]
    f.close()
    all_clusters = {}
    cluster = []
    for i in range(len(lines)):
        if i==0:
            name = lines[i][1:].split("|")
            name = "_".join(name)
        elif i==len(lines)-1:
            cluster.append(lines[i])
            if "" not in cluster: #there might be clusters without sequences
                all_clusters[name]=cluster
                with open(f"{outclusterpath}/{name}.fasta", "w") as f:
                    for c in cluster:
                        f.write(c+"\n")
                    f.close()
        else:
            if lines[i]==lines[i+1]:
                if "" not in cluster: #there might be clusters without sequences
                    all_clusters[name] = cluster
                    with open(f"{outclusterpath}/{name}.fasta", "w") as f:
                        for c in cluster:
                            f.write(c+"\n")
                        f.close()
                name = lines[i][1:].split("|")
                name = "_".join(name)
                cluster = []
            else:
                cluster.append(lines[i])
                
    # filter only clusters with 38 genes and 1-1 
    clusters_filtered = {}
    for name in all_clusters:
        if len(all_clusters[name])/2 == len(genomeidlist):
            names_in_cluster = set()
            unwanted_aa = False
            for c in all_clusters[name]:
                if c[0]==">":
                    names_in_cluster.add(c[1:].split("|")[0])
                else:
                    # check if no J|U|O aminoacids
                    if search_unwanted_aminoacids(c): #found unwanted 
                        unwanted_aa = True
            if (set(genomeidlist)==names_in_cluster) and (unwanted_aa==False):
                clusters_filtered[name] = all_clusters[name].copy()
                #save filtered
                with open(f"{filteredclusterpath}/{name}.fasta", "w") as f:
                    for c in all_clusters[name]:
                        if c[0]==">": #rename gene name to genome name
                            f.write(c.split("|")[0]+"\n")
                        else:
                            f.write(c+"\n")
                    f.close()
    # PARALOGS
    if paralogsclusterpath:
        if not os.path.exists(paralogsclusterpath):
            os.makedirs(paralogsclusterpath)
        # filter clusters allowing paralogs, but at least 3 different genomes
        clusters_paralogs = {}
        for name in all_clusters:
            if len(all_clusters[name])/2 >= 3: #at least 3 sequences
                names_in_cluster = set()
                unwanted_aa = False
                for c in all_clusters[name]:
                    if c[0]==">":
                        names_in_cluster.add(c[1:].split("|")[0])
                    else:
                        # check if no J|U|O aminoacids
                        if search_unwanted_aminoacids(c): #found unwanted 
                            unwanted_aa = True
                if (len(names_in_cluster) >= 3) and (unwanted_aa==False): #at least 3 different genomes present
                    clusters_paralogs[name] = all_clusters[name].copy() 
                    with open(f"{paralogsclusterpath}/{name}.fasta", "w") as f:
                        for c in all_clusters[name]:
                            f.write(c+"\n")
                        f.close()
    #     return all_clusters, clusters_filtered, clusters_paralogs
    # else: # WITHOUT paralogs          
    #     return all_clusters, clusters_filtered


##########################################################################################################################
##########################################################################################################################

# Clustal W - alignment
def align_clustalw(clusterpath, alignpath):
    if not os.path.exists(alignpath):
        os.makedirs(alignpath)
    # alignment with ClustalW
    for in_file in glob(f"{clusterpath}/*.fasta", recursive=True): 
        name = in_file.split("\\")[-1][:-6]
        out_file = f"{alignpath}/{name}.aln"
        clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file, outfile=out_file)
        clustalw_cline()

##########################################################################################################################
##########################################################################################################################

# Build NJ trees

def has_positive_branch_lengths(tree):
    for t in tree.get_terminals():
        if t.branch_length < 0:
            return False
    for t in tree.get_nonterminals():
        if t.branch_length < 0 :
            return False
    return True

def get_mean_support(target_tree, bootstrap_trees):
    support_tree = get_support(target_tree, bootstrap_trees)
    count = 0
    conf_sum = 0
    for t in support_tree.get_nonterminals():
        if t.confidence:
            count+=1
            conf_sum+=t.confidence
    return conf_sum/count

##########################################################################################################################

def calculate_nj_trees(alignpath, treespath, treesurecpath, bootstrappath=None, mean_support=75, bn= 50):

    # create directories
    if not os.path.exists(treespath):
        os.makedirs(treespath)
    if not os.path.exists(treesurecpath):
        os.makedirs(treesurecpath)

    if bootstrappath:
        bootstraptreepath = bootstrappath + f"/trees_{mean_support}"
        if not os.path.exists(bootstraptreepath):
            os.makedirs(bootstraptreepath)

        bootstraptreeurecpath = bootstrappath + f"/trees_urec_{mean_support}"
        if not os.path.exists(bootstraptreeurecpath):
            os.makedirs(bootstraptreeurecpath)

        #filtertrees
        filtered_trees = []

    #calculate trees
    calculator = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor()
    trees = []

    for path in glob(f"{alignpath}/*.aln", recursive=True):
        #read alignment
        aln = AlignIO.read(path, 'clustal')
        distance_matrix = calculator.get_distance(aln)
        target_tree = constructor.nj(distance_matrix)
        # save only if all branch lenghts are positive
        if has_positive_branch_lengths(target_tree) == True:
            #save name
            name = path.split("\\")[-1][:-4]
            trees.append(target_tree)
            #out file
            out_file = f"{treespath}/{name}.nwk"
            #save the tree
            Phylo.write(target_tree, out_file, "newick")
            #save with no internals
            for t in target_tree.get_nonterminals():
                t.name=None
            out_file = f"{treesurecpath}/{name}.nwk"
            Phylo.write(target_tree, out_file, "newick")
            #########################################################################
            if bootstrappath:
                #bootstrap
                aln_reps = bootstrap(aln, bn)
                bootstrap_trees = []
                for a in aln_reps:
                    distance_matrix = calculator.get_distance(a)
                    NJ_tree = constructor.nj(distance_matrix)
                    #save if has positive branch lengths
                    if has_positive_branch_lengths(NJ_tree) == True:
                        bootstrap_trees.append(NJ_tree)
                # calculate support of target tree and save if above threshold
                if get_mean_support(target_tree, bootstrap_trees)>=mean_support:
                    name = path.split("\\")[-1][:-4]
                    filtered_trees.append(target_tree)
                    #out file
                    out_file = f"{bootstraptreepath}/{name}.nwk"
                    #save the tree
                    Phylo.write(target_tree, out_file, "newick")
                    #save with no internals
                    for t in NJ_tree.get_nonterminals():
                        t.name=None
                    out_file = f"{bootstraptreeurecpath}/{name}.nwk"
                    Phylo.write(target_tree, out_file, "newick")
            #########################################################################

    if bootstrappath:
        return trees, filtered_trees, bootstraptreepath, bootstraptreeurecpath
    else:
        return trees

##########################################################################################################################

def calculate_nj_trees_paralogs(alignparapath, paralogspath):
    
    paralogstreepath = paralogspath + f"/trees"
    if not os.path.exists(paralogstreepath):
        os.makedirs(paralogstreepath)

    paralogstreeurecpath = paralogspath + f"/trees_urec"
    if not os.path.exists(paralogstreeurecpath):
        os.makedirs(paralogstreeurecpath)
    
    calculator = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor()

    para_trees = []

    for path in glob(f"{alignparapath}/*.aln", recursive=True):
        #read alignment
        aln = AlignIO.read(path, 'clustal')
        distance_matrix = calculator.get_distance(aln)
        NJ_tree = constructor.nj(distance_matrix)
        # save only if all branch lenghts are positive
        if has_positive_branch_lengths(NJ_tree) == True:
            #save name
            name = path.split("\\")[-1][:-4]
            para_trees.append(NJ_tree)
            #out file
            out_file = f"{paralogstreepath}/{name}.nwk"
            #save the tree
            Phylo.write(NJ_tree, out_file, "newick")
            #save with no internals
            for t in NJ_tree.get_nonterminals():
                t.name = None
            #save with only genome names 
            for t in NJ_tree.get_terminals():
                t.name = t.name.split("|")[0]
            out_file = f"{paralogstreeurecpath}/{name}.nwk"
            Phylo.write(NJ_tree, out_file, "newick")
            
    return para_trees, paralogstreepath, paralogstreeurecpath

##########################################################################################################################
##########################################################################################################################

# FASTUREC, supertrees

def prepare_for_fasturec(treesurecpath, treeforurecpath, result_file):
    urec_trees = []
    for path in glob(f"{treesurecpath}/*.nwk", recursive=True):
        tree = Tree(path, format=1)
        name = path.split("/")[-1][:-4]
        out_file = f"{treesurecpath}/{name}.nwk"
        tree.write(outfile=out_file, format=9)
        urec_trees.append(tree)

    with open(f"{treeforurecpath}/{result_file}.nwk", "w") as f:
        for t in urec_trees:
            t_str = t.write(format=9)
            t_str = t_str.replace('NoName', '')
            f.write(t_str[:-1]+"\n")

def get_supertrees(treeforurecpath, bootstrap=False, paralogs=False):
    
    #run 1-1 inference
    subprocess.call(f"./fasturec -G {treeforurecpath}/tree_urec.nwk -Y", shell=True)

    if bootstrap: # run bootstrap inference
        subprocess.call(f"./fasturec -G {treeforurecpath}/tree_urec_bootstrap.nwk -Y", shell=True)
    
    if paralogs: # run paralogs inference
        subprocess.call(f"./fasturec -G {treeforurecpath}/tree_urec_paralogs.nwk -Y", shell=True)


def save_optimal_supertree(fasturecpath, supertreepath, in_file, out_file):
    for file in os.listdir(fasturecpath):
        if file.startswith(f"fu.{in_file}."):
            supertrees = []
            optimal_cost = file.split(".")[-3]
            print(f"Optimal cost of supertree: {optimal_cost}")
            with open(file, "r") as f:
                for l in f.readlines():
                    line = l.strip().split(" ")
                    if line[0]==optimal_cost:
                        supertrees.append(line[1])
            print(f"We have {len(supertrees)} optimal trees")
    # save results
    for i in range(len(supertrees)):
        super_tree = Tree(supertrees[i]+";")
        super_tree.write(outfile=f"{supertreepath}/{out_file}_{i}.nwk", format=9)

##########################################################################################################################
##########################################################################################################################

def save_for_consensus(trees, consensustreepath, out_file):
    Phylo.write(trees, f"{consensustreepath}/{out_file}.nwk", "newick")

def get_consensus(bootstrap=False):

    # 1-1 inference
    subprocess.call('/usr/bin/Rscript consensus_tree.R', shell=True)
    # bootstrap inference
    if bootstrap:
        subprocess.call('/usr/bin/Rscript consensus_tree.R -b', shell=True)

##########################################################################################################################
##########################################################################################################################

def main():

    ###########################################################################################################################
    args = parse_arguments()

    MAINDIR = os.getcwd()
    DATADIR = MAINDIR+"/data"
    if args.email:
        EMAIL = args.email
    else:
        EMAIL = "m.bochenek@student.uw.edu.pl"

    GENOMELISTFILE = MAINDIR + f"/{args.n}"
    
    # number of genomes
    genomeidlist = processGenomeList(GENOMELISTFILE)

    ###########################################################################################################################
    # Download data

    # Create directory for genbank files storage
    datapath = DATADIR
    if(len(glob(datapath+"\*")) > 0):
        os.system("rm -rf "+datapath+"\*")

    gbkpath = datapath + "/gbk"
    if not os.path.exists(gbkpath):
        os.makedirs(gbkpath)
    os.chdir(gbkpath)

    # retrieve genomes
    retrieve_genomes(EMAIL, GENOMELISTFILE)

    #retrieve genes
    genespath = datapath + "/genes"
    if not os.path.exists(genespath):
        os.makedirs(genespath)

    download_genes(gbkpath, genespath+"/genes.fasta")

    ###########################################################################################################################
    # Cluster genes

    # call cluster_genes.sh
    os.chdir(MAINDIR)
    subprocess.call("./cluster_genes.sh", shell=True)

    ###########################################################################################################################
    # Reading clusters and Alignment 

    #mmseq cluster path
    clusterpath = datapath + "/clustering"
    #out cluster path
    outclusterpath = datapath + "/out_clusters"
    # filtered clusters "1-1"
    filteredclusterpath = datapath + "/filtered_clusters"

    if args.p:
        # filtered clusters with paralogs and more than 3 sequences
        paralogsclusterpath = datapath + "/paralogs_clusters"
        # filter clusters
        read_and_filter_clusters(genomeidlist, clusterpath, outclusterpath, filteredclusterpath, paralogsclusterpath)
        #Alignment 
        alignparapath = datapath + "/align_paralogs"
        # align
        align_clustalw(paralogsclusterpath, alignparapath)
    else:
        read_and_filter_clusters(genomeidlist, clusterpath, outclusterpath, filteredclusterpath)

    #Alignment - ClustalW
    alignpath = datapath + "/align"
    # align
    align_clustalw(filteredclusterpath, alignpath)

    ###########################################################################################################################
    # NJ trees

    # Build NJ trees
    treespath = datapath + "/trees"
    treesurecpath = datapath + "/trees_urec"

    # folder for bootstrapped NJ trees
    if args.b:
        bootstrappath = datapath + "/bootstrap"
        if not os.path.exists(bootstrappath):
            os.makedirs(bootstrappath)

    #save trees for consensus calculation
    consensustreepath = datapath + "/consensus"
    if not os.path.exists(consensustreepath):
        os.makedirs(consensustreepath)

    # calculate nj trees 
    if args.b:
        nj_trees, bootstrap_trees, bootstraptreepath, bootstraptreeurecpath = calculate_nj_trees(alignpath, treespath, treesurecpath, bootstrappath=bootstrappath, mean_support=args.thr, bn= args.bn)
        #save tree for consensus calculation
        save_for_consensus(bootstrap_trees, consensustreepath, out_file="trees_for_consensus_bootstrap")
    else:
        nj_trees = calculate_nj_trees(alignpath, treespath, treesurecpath)

    # save all the nj trees to one file for consensus calculation 
    save_for_consensus(nj_trees, consensustreepath, out_file="trees_for_consensus")

    # PARALOGS - NJ trees
    # folder for paralogs NJ trees
    if args.p:
        paralogspath = datapath + "/paralogs"
        if not os.path.exists(paralogspath):
            os.makedirs(paralogspath)
        # calculate NJ trees
        para_trees, paralogstreepath, paralogstreeurecpath = calculate_nj_trees_paralogs(alignparapath, paralogspath)

    ###########################################################################################################################
    # Supertrees

    # FASTUREC
    #prepare trees for fasturec
    treeforurecpath = datapath + "/tree_for_urec"
    if not os.path.exists(treeforurecpath):
        os.makedirs(treeforurecpath)

    #prepare for fasturec
    prepare_for_fasturec(treesurecpath, treeforurecpath, "tree_urec")
    if args.b:
        prepare_for_fasturec(bootstraptreeurecpath, treeforurecpath, "tree_urec_bootstrap")
    if args.p:
        prepare_for_fasturec(paralogstreeurecpath, treeforurecpath, "tree_urec_paralogs")

    # get supertrees 
    fasturecpath = MAINDIR 
    os.chdir(fasturecpath)

    # get supertress
    get_supertrees(treeforurecpath, bootstrap=args.b, paralogs=args.p)

    # save optimal super trees
    supertreepath = datapath + "/supertree"
    if not os.path.exists(supertreepath):
        os.makedirs(supertreepath)

    save_optimal_supertree(fasturecpath, supertreepath, in_file="tree_urec.nwk", out_file="super_tree")
    if args.b:
        save_optimal_supertree(fasturecpath, supertreepath, in_file="tree_urec_bootstrap.nwk", out_file="super_tree_bootstrap")
    if args.p:
        save_optimal_supertree(fasturecpath, supertreepath, in_file="tree_urec_paralogs.nwk", out_file="super_tree_paralogs")

    ###########################################################################################################################

    # Consensus tree - withount and with bootstraping
    get_consensus(bootstrap=args.b)

##########################################################################################################################
##########################################################################################################################

if __name__ == "__main__":
    main()
