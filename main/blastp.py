#!/usr/bin/python
# coding: utf-8

import os
import sys
import time
from Bio import SeqIO


def main() :
    start_time = time.time()
    name = sys.argv[1:][0]
    genes = list()
    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            genes.append(gene)

    for gene in genes :
        print(gene)
        # Execute shell commands in Python
        myCmd = "blastp -db nr -remote -query ./%s -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue length pident' -out ./%s.txt" % (name, gene)
        os.system(myCmd)

    elapsed_time = (time.time() - start_time)
    print("Finished------------------")
    print("Running the Blast process: %.4fs" %(elapsed_time))


if __name__ == "__main__" :
    main()