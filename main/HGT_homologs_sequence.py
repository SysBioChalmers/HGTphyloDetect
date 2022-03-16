#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
# This file is to get 300 homolog accession id and sequences for the following multiple sequence alignment.


import re
import sys
import os
from Bio import SeqIO
from Bio import Entrez
from ete3 import NCBITaxa
from subprocess import Popen, PIPE


LINNAEUS_FILTER = ["species","genus","family","order","class","subphylum","phylum","kingdom","superkingdom"]

def parse_NCBI(filename):
    with open(filename, "r") as infile :
        lines = infile.readlines()

    accession_number = list()
    accession_similarity = dict()

    for line in lines :
        accession = line.strip("\n").split("\t")[2]
        similarity = line.strip("\n").split("\t")[3]
        accession_number.append(accession)
        accession_similarity[accession] = float(similarity)

    # print(len(accession_number))
    return accession_number, accession_similarity

def get_refSeq(gene) :
    # get the protein sequence accoding to protein sequence id
    with open("./fasta/%s.fasta" % gene, "r") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            proteinSeq[record.id] = str(record.seq)

    return proteinSeq[gene]

accession_taxid = dict()
with open("acc2taxid_final.txt","r") as infile :
    lines = infile.readlines()
for line in lines :
    accession = line.strip('\n').split(',')[0]
    taxid = line.strip('\n').split(',')[1]
    accession_taxid[accession] = taxid

def getTaxid(accession) :
    taxid = accession_taxid[accession]
    return taxid

def getSeq(accession):
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

    seq = record.seq

    print(seq)

    return seq

def main() : 
    with open("saccharomyces_cerevisiae_S288c_HGT.tsv", "r") as infile :
        lines = infile.readlines()[1:]

    genes = [line.strip().split('\t')[0] for line in lines]
    # print(genes)
    # print(len(genes))  # 9

    i = 0
    for gene in genes :
        i += 1
        print('This is', str(i) + '--------------------------------------------------------')
        print(gene)
        accession_number, accession_similarity = parse_NCBI("./blastp_files/%s.txt" % gene)
        # print(len(accession_number)) # 300

        ncbi = NCBITaxa()
        id_seq = dict()
        for accession in accession_number :
            # print(accession)
            # print(accession_similarity[accession])
            try :
                seq = str(getSeq(accession))
            except :
                continue

            taxid = getTaxid(accession)
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            # print(ranks2lineage)

            try :
                taxid2name = ncbi.get_taxid_translator([ranks2lineage['phylum'], ranks2lineage['species']])
                if taxid2name[ranks2lineage["phylum"]]  != "Ascomycota":
                    accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@" + taxid2name[ranks2lineage["phylum"]]

                # if taxid2name[ranks2lineage["superkingdom"]]  == "Bacteria":
                #     print(taxid2name[ranks2lineage["phylum"]])
                #     print(taxid2name[ranks2lineage["superkingdom"]])

                    print(accession)
                    id_seq[accession] = seq
                    continue
            except :
                pass

            try :
                taxid2name = ncbi.get_taxid_translator([ranks2lineage['phylum'], ranks2lineage['subphylum'], ranks2lineage['species']])
                # print(taxid2name[ranks2lineage["subphylum"]])

                # The Eurotiomycetes are a class of ascomycetes within the subphylum Pezizomycotina.
                # Proteobacteria is a major phylum of Gram-negative bacteria.
                if accession_similarity[accession] < 80 :  # Sequences with more than 80% similarity were eliminated.
                    if taxid2name[ranks2lineage["subphylum"]]  == "Saccharomycotina" :
                        accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@Saccharomycotina"
                    if taxid2name[ranks2lineage["subphylum"]]  != "Saccharomycotina" and taxid2name[ranks2lineage["phylum"]]  == "Ascomycota":
                        accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@other_Ascomycota"
                    if taxid2name[ranks2lineage["phylum"]]  != "Ascomycota":
                        accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@" + taxid2name[ranks2lineage["phylum"]]
                    # if taxid2name[ranks2lineage["kingdom"]]  != "Fungi"  :
                    #     accession = accession + "-" + taxid2name[ranks2lineage["species"]].replace(" ", "_") + "@Other"
                    print(accession)
                    id_seq[accession] = seq
                else :
                    continue
            except :
                continue

        gene_seq = get_refSeq(gene)
        gene_query = "QUERY_" + 'Saccharomyces_cerevisiae_' + gene
        print(gene_query)
        print(gene_seq)
        id_seq[gene_query] = gene_seq
        print(len(id_seq))

        os.system("mkdir ./fasta_homologs/%s" % gene)
        with open("./fasta_homologs/%s/%s.fasta" % (gene,gene), "w") as outfile :
            for accession,seq in id_seq.items() :
                outfile.write(">"+accession+"\n")
                outfile.write(seq+"\n")


if __name__== "__main__":
    main()
