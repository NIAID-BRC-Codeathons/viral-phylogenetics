#!/bin/bash

mkdir -p proteome/{species,filtered}

# generate full FASTA file to start with
mmseqs convert2fasta pdb_virus_foldseek/pdb_virus proteome/pdb_virus.fa

# use tax info to split by species
awk 'NR==FNR{f[">"$1]=$2; next} $1 in f{x=f[$1]; print > "proteome/species/"x".fa"; getline; print > "proteome/species/"x".fa"}' pdb_virus.taxinfo proteome/pdb_virus.fa

# get species with at least 50 sequences
wc -l proteome/species/* | awk '$1>99 && !/total/ && !/0\.fa$/ {print $2}' | cut -d/ -f3 | while read F; do cp -v proteome/species/$F proteome/filtered/$F; done
