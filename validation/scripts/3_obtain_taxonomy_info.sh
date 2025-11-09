#!/bin/bash

# To obtain taxonomy information, we need to convert structure DBs into result DBs.
# For this, we tweak mmseqs/foldseek to self align their sequences and produce result DBs.

mmseqs tsv2db <(awk '{print $1"\t"$1}' pdb_virus_foldseek/pdb_virus.index) pdb_virus_foldseek/self_res --output-dbtype 7

foldseek align pdb_virus_foldseek/pdb_virus pdb_virus_foldseek/pdb_virus pdb_virus_foldseek/self_res pdb_virus_foldseek/self_aln

# Now we can use the alignment DB to obtain taxonomy information.
mmseqs convertalis pdb_virus_foldseek/pdb_virus pdb_virus_foldseek/pdb_virus pdb_virus_foldseek/self_aln pdb_virus.taxinfo --format-output query,taxid,taxname,taxlineage
