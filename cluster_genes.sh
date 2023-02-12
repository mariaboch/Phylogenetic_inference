#!/bin/bash
cd data
mkdir clustering
mmseqs createdb genes/genes.fasta clustering/DB
mmseqs cluster clustering/DB clustering/DB_clu tmp
mmseqs createtsv clustering/DB clustering/DB clustering/DB_clu clustering/DB_clu.tsv
mmseqs createseqfiledb clustering/DB clustering/DB_clu clustering/DB_clu_seq
mmseqs result2flat clustering/DB clustering/DB clustering/DB_clu_seq clustering/DB_clu_seq.fasta