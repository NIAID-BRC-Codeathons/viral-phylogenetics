#!/bin/bash

mkdir -p unicore/{core-prostt5,core-pdb,model}
TMP=$(mktemp -d)

# run unicore with given proteomes, using default ProstT5 model
unicore easy-core proteome/filtered unicore/core-prostt5/ unicore/model/ $TMP

# run unicore with PDB structures
unicore easy-core proteome/filtered unicore/core-pdb/ unicore/model/ $TMP/unicore --custom-lookup pdb_virus_foldseek/pdb_virus
