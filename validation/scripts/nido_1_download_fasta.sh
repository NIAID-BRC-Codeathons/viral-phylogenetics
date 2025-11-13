#!/bin/bash

mkdir -p nido_fasta
awk -F',' 'NR==FNR{f[$1]=1;next} $1 in f{print "wget -q -O- https://rest.uniprot.org/uniprotkb/"$3".fasta >> nido_fasta/"$1".fa"}' <(cut -d, -f1 nidovirales_ref-proteomes_proteins.csv | tail -n+2 | uniq -c | awk '$1>9{print $2}') nidovirales_ref-proteomes_proteins.csv | bash
