#!/bin/bash

# rename fasta files in given directory into clusters compatible with ASTRAL
if [ "$#" -ne 1 ]; then
	echo "Usage: $0 <directory>"
	exit 1
fi

DIR=$1

mkdir -p $DIR/renamed/{aa,3di}

CUR="1"
ls $DIR/tree/fasta | while read C; do
    awk -v CUR=$CUR '/^>/{print $1"_"CUR; next} {print}' $DIR/tree/fasta/$C/aa.fasta > $DIR/renamed/aa/cluster_${CUR}_aa.fasta
	awk -v CUR=$CUR '/^>/{print $1"_"CUR; next} {print}' $DIR/tree/fasta/$C/3di.fasta > $DIR/renamed/3di/cluster_${CUR}_3di.fasta
	let CUR=CUR+1
done

