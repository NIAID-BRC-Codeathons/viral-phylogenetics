#!/bin/bash

# pdb_virus.id contains a list of PDB IDs of viral proteins
awk 'NR==FNR{f[tolower($1)]=1; next} substr($1,1,4) in f{print}' pdb_virus.id <(cat pdb_foldseek/pdb_h | tr -d '\0') > pdb_virus.header

# convert headers into compatible names for foldseek
awk '{sub(/\.cif\.gz/,"",$1); print $1}' pdb_virus.header > pdb_virus.lookup

# subset the original pdb database to create a viral protein pdb database
scripts/subdb foldseek pdb_virus.lookup pdb_foldseek/pdb pdb_virus_foldseek/pdb_virus
