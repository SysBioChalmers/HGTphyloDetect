#!/bin/bash


echo 'Begin to run!!!'

Rscript midpoint_tree.R YFR055W
perl create_iTOL_config.pl ../fasta/YFR055W/YFR055W_midpoint.tree

Rscript midpoint_tree.R YJL218W
perl create_iTOL_config.pl ../fasta/YJL218W/YJL218W_midpoint.tree

Rscript midpoint_tree.R YOL164W
perl create_iTOL_config.pl ../fasta/YOL164W/YOL164W_midpoint.tree

Rscript midpoint_tree.R YNR057C
perl create_iTOL_config.pl ../fasta/YNR057C/YNR057C_midpoint.tree

Rscript midpoint_tree.R YNR058W
perl create_iTOL_config.pl ../fasta/YNR058W/YNR058W_midpoint.tree

Rscript midpoint_tree.R YDR540C
perl create_iTOL_config.pl ../fasta/YDR540C/YDR540C_midpoint.tree

Rscript midpoint_tree.R YKL216W
perl create_iTOL_config.pl ../fasta/YKL216W/YKL216W_midpoint.tree

Rscript midpoint_tree.R YJL217W
perl create_iTOL_config.pl ../fasta/YJL217W/YJL217W_midpoint.tree

Rscript midpoint_tree.R YMR090W
perl create_iTOL_config.pl ../fasta/YMR090W/YMR090W_midpoint.tree

echo 'Yep, finish!!!'

