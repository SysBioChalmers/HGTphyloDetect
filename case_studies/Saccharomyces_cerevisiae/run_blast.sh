#!/usr/bin/env bash

diamond blastp --db /Users/leyu/Documents/Le/Databases/nr -q ./saccharomyces_cerevisiae_S288c.fasta -o ./blastp/blastp_all.txt -k 300 -f 6 qseqid qlen sseqid pident length evalue bitscore 

