#!/bin/bash

mkdir -p peplo_fasta
awk -F',' 'NR==FNR{f[$1]=1;next} $1 in f{print "wget -O- https://rest.uniprot.org/uniprotkb/"$3".fasta >> peplo_fasta/"$1".fa"}' <(cut -d, -f1 data/duplodnaviria-heunggongvirae-peploviricota-proteomes_proteins.csv | tail -n+2 | uniq -c | awk '$1>9{print $2}') data/duplodnaviria-heunggongvirae-peploviricota-proteomes_proteins.csv | bash
