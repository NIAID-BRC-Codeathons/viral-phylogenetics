#!/bin/bash

TMP=$(mktemp -d)

mkdir -p pdb_foldseek
foldseek databases PDB pdb_foldseek/pdb $TMP
rm -rf $TMP

mkdir -p pdb_virus_foldseek
