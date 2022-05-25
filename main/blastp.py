#!/usr/bin/python
# coding: utf-8

### Usage: python blastp.py input.fasta

import os
import sys
import time
from Bio import SeqIO


def main() :
    start_time = time.time()
    name = sys.argv[1:][0]
    genes = list()
    geneSeq = dict()
    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            sequence = str(record.seq)
            geneSeq[gene] = sequence
            genes.append(gene)

    for gene in genes :
        print('This is gene:', gene)
        with open('./%s.fasta' % gene, 'w') as outfile :
            outfile.write('>'+gene+'\n')
            outfile.write(geneSeq[gene]+'\n')
        if not os.path.exists("./blastp_files/") :
            os.system('mkdir ./blastp_files')
        # Execute shell command line in Python using the os.system function
        myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" % (gene, gene)
        os.system(myCmd)
        os.remove('./%s.fasta' % gene)

    elapsed_time = (time.time() - start_time)
    print("Finished------------------")
    print("Running the Blast process: %.4fs" %(elapsed_time))


if __name__ == "__main__" :
    main()