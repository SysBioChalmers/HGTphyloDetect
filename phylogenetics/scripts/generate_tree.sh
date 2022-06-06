#!/bin/bash

var=$1
gene=${var%.fas*}
echo 'This is gene '$gene' ' 
echo 'Begin to run!!!'

cd  ../

python scripts/HGT_homologs_sequence.py input/${gene}.fasta

mafft --thread 6 --auto ./input/${gene}_homologs.fasta > ./intermediate/${gene}_aln.fasta

trimal -in ./intermediate/${gene}_aln.fasta -out ./intermediate/${gene}_aln_trimmed.fasta -automated1

iqtree -nt 6 -st AA -s ./intermediate/${gene}_aln_trimmed.fasta -m TEST -mrate G4 -keep-ident -bb 1000 -pre ./intermediate/${gene}

Rscript scripts/midpoint_tree.R ${gene}

perl scripts/create_iTOL_config.pl ./intermediate/${gene}_midpoint.tree

echo 'Yep, finish!!!'

