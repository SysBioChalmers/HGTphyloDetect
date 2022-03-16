#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Step : BLAST hits were parsed to retrieve associated taxonomic information 
# 		using the NCBI's taxonomy database.

import re
from ete3 import NCBITaxa
import sys
import os
import math
import csv
from Bio import SeqIO
from Bio import Entrez
from subprocess import Popen, PIPE


PATH="./"
# LINNAEUS_FILTER = ["species","genus","family","order","class","subphylum","phylum","kingdom","superkingdom"]

LINNAEUS_FILTER = ["subphylum","kingdom","superkingdom"]

def parse_NCBI(gene):
    with open("../case/saccharomyces_cerevisiae/blastp/%s.txt" % gene, "r") as filename :
        data = filename.readlines()
        for line in data :
            if line.strip("\n").endswith("found") :
                index = data.index(line)

        blast_results = data[index+1:-1] 
        # print(len(blast_results))  # 250 200

    accession_number = list()
    evalue = dict()

    for blast in blast_results :
        accession = blast.strip("\n").split("\t")[1]
        accession_number.append(accession)
        evalue[accession] = float(blast.strip('\n').split("\t")[-3])

    # print(accession_number)
    # print(evalue)
    return accession_number, evalue

def getTaxid(accession):
    # Retrieving data in the GenBank using only the GenBank code accession in biopython
    # https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    # https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    Entrez.email = "leyu@example.org"

    # https://www.biostars.org/p/304175/
    # get tax id using only the GenBank code accession in biopython
    # handle = Entrez.efetch(db='protein', id="NP_012706.1", rettype='gb')
    handle = Entrez.efetch(db='protein', id=accession, rettype='gb')
    record = SeqIO.read(handle,'genbank')
    # print(record.features[0].qualifiers)
    if record.features[0].qualifiers['db_xref'][0].split(":")[0] == 'taxon':
        taxid = record.features[0].qualifiers['db_xref'][0].split(":")[1] # the type is a string
        organism = record.features[0].qualifiers['organism'][0]
    # seq = record.seq
    # print(taxid,organism)

    return taxid,organism
    # return seq


if __name__== "__main__":
    genes = list()
    HGT = list()
    for filename in os.listdir("../case/saccharomyces_cerevisiae/fasta") :
        # print("This is %d------------------" %(i))
        gene = filename[:-6]
        genes.append(gene)

    n=0
    for gene in genes :
        n += 1
        print("This is %d------------------" %(n))
        print(gene)
        accession_number,evalue = parse_NCBI(gene)

        i=0
        k=0
        ncbi = NCBITaxa()
        outgroup = list()
        ingroup = list()
        outgroup_dict = dict()
        ingroup_dict = dict()
        e_minus = 1e-200

        for accession in accession_number[:200] :
            i+=1
            # taxid, organism = getTaxid(accession)
            try :
                taxid,organism = getTaxid(accession)
                lineage = ncbi.get_lineage(taxid)
                lineage2ranks = ncbi.get_rank(lineage)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            except :
                k+=1
                # print('------------------------------',k)
                continue
            # taxonomy = os.system('ete3 ncbiquery --search %d --info' % int(taxid))
            # lineage = ncbi.get_lineage(taxid)
            # lineage2ranks = ncbi.get_rank(lineage)
            # ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())


            print(i)
            print(ranks2lineage)

            taxonomy_alignment = ranks2lineage
            LINNAEUS_FILTER = ["subphylum","kingdom","superkingdom"]

            # {'no rank': 2608891, 'superkingdom': 2, 'phylum': 1239, 'genus': 1350, 'family': 81852, 'class': 91061, 
            # 'order': 186826, 'species': 1857215}
            # Traceback (most recent call last):
            #   File "xxx.py", line 226, in <module>
            #     if taxonomy_alignment['kingdom'] != 4751 :  # 'Fungi': 4751   'Saccharomycotina': 147537
            # KeyError: 'kingdom'
            try :

                if taxonomy_alignment['kingdom'] != 4751 :  # 'Fungi': 4751   'Saccharomycotina': 147537
                    outgroup.append(accession)              # 'Eukaryota': 2759

                if taxonomy_alignment['kingdom'] == 4751 and taxonomy_alignment['subphylum'] != 147537 :
                    ingroup.append(accession)

            except :

                outgroup.append(accession)

        for accession_id in outgroup :
            outgroup_dict[accession_id] = evalue[accession_id]

        for accession_id in ingroup :
            ingroup_dict[accession_id] = evalue[accession_id]

        if outgroup_dict :
            min_outgroup_key = min(outgroup_dict,key=outgroup_dict.get)
            min_outgroup_evalue = outgroup_dict[min_outgroup_key]
        else :
            min_outgroup_evalue = 1

        if ingroup_dict :
            min_ingroup_key = min(ingroup_dict,key=ingroup_dict.get)
            min_ingroup_evalue = ingroup_dict[min_ingroup_key]
        else :
            min_ingroup_evalue = 1

        alienindex = math.log(min_ingroup_evalue+e_minus, math.e)-math.log(min_outgroup_evalue+e_minus, math.e)

        print(outgroup_dict)
        print(ingroup_dict)
        print(min_outgroup_evalue)
        print(min_ingroup_evalue)
        print(alienindex)

        if alienindex > 45 :
            print('This is a HGT event')
            print('Accession_id: %s' % min_outgroup_key)
            taxid, organism = getTaxid(min_outgroup_key)
            print('This HGT is from: %s' % organism)
            hgt = 'Yes'
        else:
            print('This is not a HGT')
            hgt = 'No'

        HGT.append([gene,alienindex,hgt])

        print(HGT)

    outfile = open("../case/saccharomyces_cerevisiae/results_HGT.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    tsv_writer.writerow(["gene id", "alien index", "HGT or not"])

    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)

    outfile.close()


