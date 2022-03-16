#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is to identify the HGT events from fungi to yeast

import re
from ete3 import NCBITaxa
import sys
import os
import math
import csv
from Bio import SeqIO
from Bio import Entrez
from subprocess import Popen, PIPE


def parse_NCBI(gene):
    with open("./blastp_files/%s.txt" % gene, "r") as filename :
        lines = filename.readlines()

    accession_number = list()
    accession_bitscore = dict()
    for line in lines :
        accession = line.strip("\n").split("\t")[2]
        accession_number.append(accession)
        bitscore = line.strip("\n").split("\t")[-1]
        accession_bitscore[accession] = float(bitscore)

    return accession_number, accession_bitscore

def getTaxid2(accession):
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

def main() :
    genes = list()
    HGT = list()

    outfile = open("./saccharomyces_cerevisiae_S288c_HGT_eukaryote2eukaryote.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    # column = ['Gene', 'Alien index', 'E value', 'Donor id', 'Donor phylum', 'Donor genus']
    column = ['Gene', 'Bitscore', 'Outg_pct', 'HGT index', 'Donor taxonomy']
    tsv_writer.writerow(column)

    for filename in os.listdir("./blastp_files") :
        # print("This is %d------------------" %(i))
        gene = filename[:-4]
        genes.append(gene)

    n=0
    k=0
    for gene in genes :
        n += 1
        print("This is %d------------------" %(n))
        # print(gene)
        accession_number, accession_bitscore = parse_NCBI(gene)

        recipient_yeast = list()
        other_fungi_accession = list()
        recipient_yeast_accession_bitscore = dict()
        other_fungi_accession_bitscore = dict()
        recipient_species = list()
        other_fungi_species = list()

        # i = 0
        ncbi = NCBITaxa()
        for accession in accession_number :
            try :
                taxid = getTaxid(accession)
                # taxid = getTaxid2(accession)
                # print(taxid)
            except :
                continue
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            # print(ranks2lineage)

            taxid2name = ncbi.get_taxid_translator(lineage)
            taxonomy_alignment = ranks2lineage

            try :
                if taxonomy_alignment['subphylum'] == 147537 :
                # if taxonomy_alignment['genus'] == 4930 :
                    recipient_yeast.append(accession)
                    recipient_species.append(taxonomy_alignment['species'])

                # if taxonomy_alignment['superkingdom'] == 2759 and taxonomy_alignment['subphylum'] != 147537 :
                if taxonomy_alignment['kingdom'] == 4751 and taxonomy_alignment['subphylum'] != 147537 :
                # if taxonomy_alignment['kingdom'] == 4751 and taxonomy_alignment['phylum'] != 4890 :  # 'Ascomycota' phylum
                    other_fungi_accession.append(accession)
                    other_fungi_species.append(taxonomy_alignment['species'])
            except :
                continue

        for accession_id in recipient_yeast :
            recipient_yeast_accession_bitscore[accession_id] = accession_bitscore[accession_id]

        for accession_id in other_fungi_accession :
            other_fungi_accession_bitscore[accession_id] = accession_bitscore[accession_id]

        if recipient_yeast_accession_bitscore :
            max_recipient_yeast_accession_key = max(recipient_yeast_accession_bitscore,key=recipient_yeast_accession_bitscore.get)
            max_recipient_yeast_bitscore = recipient_yeast_accession_bitscore[max_recipient_yeast_accession_key]

        if other_fungi_accession_bitscore :
            max_other_fungi_accession_key = max(other_fungi_accession_bitscore,key=other_fungi_accession_bitscore.get)
            max_other_fungi_bitscore = other_fungi_accession_bitscore[max_other_fungi_accession_key]
            if max_other_fungi_accession_key :
                max_taxid = getTaxid(max_other_fungi_accession_key)
                max_lineage = ncbi.get_lineage(max_taxid)
                max_lineage2ranks = ncbi.get_rank(max_lineage)
                max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
                try :
                    max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage['kingdom'], max_ranks2lineage['phylum'], max_ranks2lineage['subphylum'], max_ranks2lineage['species']])
                except :
                    continue

                print(gene)
                print(max_recipient_yeast_bitscore)
                print(max_other_fungi_bitscore)
                print(max_taxid2name)

                if recipient_species :
                    recipient_species_number = len(set(recipient_species))
                if other_fungi_species :
                    other_fungi_species_number = len(set(other_fungi_species))

                HGT_index = format(max_other_fungi_bitscore/max_recipient_yeast_bitscore, '.4f')
                Outg_pct = format(other_fungi_species_number/(other_fungi_species_number+recipient_species_number), '.4f')

                print(HGT_index)
                print(recipient_species_number)
                print(other_fungi_species_number)
                print(Outg_pct) 
                # if max_other_fungi_bitscore>100 and float(HGT_index)>=0.5 and float(Outg_pct)>=0.9 :
                if max_other_fungi_bitscore>100 and float(HGT_index)>=0.5 and float(Outg_pct)>=0.8 :
                    k += 1 
                    print("This is a potential HGT event from fungi to yeast!!!")
                    print('This is HGT',k)

                    taxonomy = max_taxid2name[max_ranks2lineage['kingdom']] + '/' + max_taxid2name[max_ranks2lineage['subphylum']]
                    tsv_writer.writerow([gene, max_other_fungi_bitscore, Outg_pct, HGT_index, taxonomy])

    outfile.close()


if __name__== "__main__":
    main()
