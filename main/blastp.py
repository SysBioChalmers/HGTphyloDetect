#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import time


def main() :
    start_time = time.time()
    i = 0
    for filename in os.listdir("../case/saccharomyces_cerevisiae/fasta") :  # Saccharomyces cerevisiae S288C
        # print(filename)
        i += 1
        print("This is gene %d------------------" %(i))
        gene = filename[:-6]
        print(gene)
        # Execute shell commands in Python
        # myCmd = "blastp -db nr -query ../case/saccharomyces_cerevisiae/fasta/%s.fasta -max_target_seqs 250 -num_threads 8 -task 'blastp-fast' -outfmt '7 qacc sacc evalue length pident' -out ../case/saccharomyces_cerevisiae/blastp/%s.txt" %(gene, gene)
        myCmd = "blastp -db nr -remote -query ../case/saccharomyces_cerevisiae/fasta/%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue length pident' -out ../case/saccharomyces_cerevisiae/blastp/%s.txt" %(gene, gene)

        os.system(myCmd)

    elapsed_time = (time.time() - start_time)
    print("Finished------------------")
    print("Running these gene targets takes: %.4fs" %(elapsed_time))


if __name__ == "__main__" :
    main()