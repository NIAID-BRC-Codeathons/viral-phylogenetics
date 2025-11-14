#!/bin/bash

# This script is a wrapper for the unicore command line tool, specifically designed to utilize precomputed 3Di sequences.
# The final output will be compatible with the subsequent analysis using foldtree.
#
# Input  : A directory containing two folders (aa, 3di), each contains set of FASTA files of proteomes with matching names.
#          Each proteome should have two files:
#		   - aa/<name>.fa  (amino acid sequences)
#		   - 3di/<name>.fa (3Di-represented structures)
# Output : A directory of structural core gene clusters.
#          Each cluster will be represented as a directory containing:
#          - A FASTA file with sequences from the core genes. (<cluster_id>.aa.fa)
#          - A FASTA file with 3Di-represented structures. (<cluster_id>.3di.fa)
#          - A Foldseek database built with above files. (<cluster_id>_db*)
#
# Usage  : ./unicore_wrapper.sh <input_dir> <output_dir> <model_dir> <tmp_dir> <unicore_params>
# e.g.   : ./unicore_wrapper.sh proteomes output model tmp "--gpu --max-len 2000 --core-threshold 50"

# Check for correct number of arguments
if [ "$#" -lt 5 ]; then
	echo "Usage: $0 <input_dir> <output_dir> <model_dir> <tmp_dir> <unicore_params>"
	exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
MODEL_DIR=$3
TMP_DIR=$4
UNICORE_PARAMS=${@:5}

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check unicore and foldseek dependency
if ! command -v unicore &> /dev/null; then
	echo "Error: unicore is not installed or not in PATH."
	exit 1
fi
if ! command -v foldseek &> /dev/null; then
	echo "Error: foldseek is not installed or not in PATH."
	exit 1
fi

# Convert sequences into Foldseek database
mkdir -p "$OUTPUT_DIR"/lookup
MMSEQS_FORCE_MERGE foldseek base:createdb "$INPUT_DIR"/aa "$OUTPUT_DIR"/lookup/db
MMSEQS_FORCE_MERGE foldseek base:createdb "$INPUT_DIR"/3di "$OUTPUT_DIR"/lookup/db_ss
foldseek base:rmdb "$OUTPUT_DIR"/lookup/db_ss_h

# Run unicore with provided parameters
unicore easy-core "$INPUT_DIR" "$OUTPUT_DIR" "$MODEL_DIR" "$TMP_DIR" $UNICORE_PARAMS --custom-lookup "$OUTPUT_DIR"/lookup/db

# Post-process each cluster to create Foldseek databases
CUR=1
mkdir -p $OUTPUT_DIR/core_genes
ls $OUTPUT_DIR/tree/fasta | while read C; do
	mkdir -p $OUTPUT_DIR/core_genes/cluster_$CUR
	awk -v CUR=$CUR '/^>/{print $1"_"CUR; next} {print}' $OUTPUT_DIR/tree/fasta/$C/aa.fasta > $OUTPUT_DIR/core_genes/cluster_$CUR/cluster_$CUR.aa.fa
	awk -v CUR=$CUR '/^>/{print $1"_"CUR; next} {print}' $OUTPUT_DIR/tree/fasta/$C/3di.fasta > $OUTPUT_DIR/core_genes/cluster_$CUR/cluster_$CUR.3di.fa
	MMSEQS_FORCE_MERGE=1 foldseek base:createdb $OUTPUT_DIR/core_genes/cluster_$CUR/cluster_$CUR.aa.fa $OUTPUT_DIR/core_genes/cluster_$CUR/db
	MMSEQS_FORCE_MERGE=1 foldseek base:createdb $OUTPUT_DIR/core_genes/cluster_$CUR/cluster_$CUR.3di.fa $OUTPUT_DIR/core_genes/cluster_$CUR/db_ss
	foldseek base:rmdb $OUTPUT_DIR/core_genes/cluster_$CUR/db_ss_h
	CUR=$((CUR + 1))
done
