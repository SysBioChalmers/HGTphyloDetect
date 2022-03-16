#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is to write shell scripts to generate the midpoint trees and then generate font and color


with open("../saccharomyces_cerevisiae_S288c_HGT.tsv", "r") as infile :
    lines = infile.readlines()[1:]

genes = [line.strip().split('\t')[0] for line in lines]
print(len(genes))  # 9
print(genes[:3])

file = open("generate_all_trees.sh", "w")

file.write("#!/bin/bash\n")
file.write("\n")
file.write("\n")

file.write("echo 'Begin to run!!!'\n")
file.write("\n")

for gene in genes :
    file.write("Rscript midpoint_tree.R %s\n" % gene)
    file.write("perl create_iTOL_config.pl ../fasta/%s/%s_midpoint.tree\n" % (gene, gene))
    file.write("\n")


file.write("echo 'Yep, finish!!!'\n")
file.write("\n")

file.close()













