#!/bin/bash


echo 'Begin to run!!!'

mafft --thread 6 --auto ./input/YOL164W_homologs.fasta > ./intermediate/YOL164W_aln.fasta

trimal -in ./intermediate/YOL164W_aln.fasta -out ./intermediate/YOL164W_aln_trimmed.fasta -automated1

iqtree -nt 6 -st AA -s ./intermediate/YOL164W_aln_trimmed.fasta -m TEST -mrate G4 -keep-ident -bb 1000 -pre ./intermediate/YOL164W

Rscript scripts/midpoint_tree.R YOL164W

perl scripts/create_iTOL_config.pl ./intermediate/YOL164W_midpoint.tree

echo 'Yep, finish!!!'

