#!/bin/bash

# downloads FASTA file from UniProt with given identifiers

mkdir -p fasta tmp
IDS=$1

cat $IDS | while read ID; do
    datasets download virus genome accession $ID --include protein --filename tmp/ncbi_dataset_$ID.zip
	unzip -o tmp/ncbi_dataset_$ID.zip -d tmp/ncbi_dataset_$ID
	cp -v tmp/ncbi_dataset_$ID/ncbi_dataset/data/protein.faa fasta/$ID.fasta
done


