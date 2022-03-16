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


# ncbi gi dump files location
PATH="./"

# LINNAEUS_FILTER = ["species","genus","family","order","class","subphylum","phylum","kingdom","superkingdom"]

LINNAEUS_FILTER = ["subphylum","kingdom","superkingdom"]

def parse_NCBI(gene):
    # with open("./blastp/Metschnikowia_borealis@Seq_758.txt", "r") as filename :
    with open("./blastp_files/%s.txt" % gene, "r") as filename :
        lines = filename.readlines()

    accession_number = list()
    evalue = dict()

    for line in lines :
        accession = line.strip("\n").split("\t")[2]
        accession_number.append(accession)
        bitscore = line.strip("\n").split("\t")[-1]
        evalue[accession] = float(line.strip("\n").split("\t")[-2])

    return accession_number, evalue

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

    outfile = open("./candida_versatilis_HGT.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    # column = ['Gene', 'Alien index', 'E value', 'Donor id', 'Donor phylum', 'Donor genus']
    column = ['Gene', 'Alien index', 'Outg_pct', 'E value', 'Donor id', 'Donor taxonomy']
    tsv_writer.writerow(column)

    for filename in os.listdir("./blastp_files") :
        # print("This is %d------------------" %(i))
        gene = filename[:-4]
        genes.append(gene)

    n=0
    for gene in genes :
        n += 1
        print("This is %d------------------" %(n))
        print(gene)
        accession_number,evalue = parse_NCBI(gene)

    # code1
        i=0
        k=0
        ncbi = NCBITaxa()
        outgroup = list()
        ingroup = list()
        outgroup_species = list()
        ingroup_species = list()
        outgroup_dict = dict()
        ingroup_dict = dict()
        e_minus = 1e-200

        for accession in accession_number :
            i+=1
            # taxid, organism = getTaxid(accession)
            try :
                taxid = getTaxid(accession)
            except :
                k+=1
                # print('------------------------------',k)
                continue
            # taxonomy = os.system('ete3 ncbiquery --search %d --info' % int(taxid))
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())


            # print(i)
            # print(ranks2lineage)

            taxonomy_alignment = ranks2lineage
            LINNAEUS_FILTER = ["subphylum","kingdom","superkingdom"]

            # {'no rank': 2608891, 'superkingdom': 2, 'phylum': 1239, 'genus': 1350, 'family': 81852, 'class': 91061, 
            # 'order': 186826, 'species': 1857215}
            # Traceback (most recent call last):
            #   File "workflow_batch_ete3.py", line 226, in <module>
            #     if taxonomy_alignment['kingdom'] != 4751 :  # 'Fungi': 4751   'Saccharomycotina': 147537
            # KeyError: 'kingdom'
            try :

                if taxonomy_alignment['kingdom'] != 4751 :  # 'Fungi': 4751   'Saccharomycotina': 147537
                    outgroup.append(accession)              # 'Eukaryota': 2759
                    outgroup_species.append(taxonomy_alignment['species'])

                if taxonomy_alignment['kingdom'] == 4751 and taxonomy_alignment['subphylum'] != 147537 :
                    ingroup.append(accession)
                    ingroup_species.append(taxonomy_alignment['species'])

            except :

                outgroup.append(accession)
                try :
                    outgroup_species.append(taxonomy_alignment['species'])
                except :
                    continue

        print('outgroup species:', len(set(outgroup_species)))
        print('ingroup species:', len(set(ingroup_species)))

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

        alienindex = format(math.log(min_ingroup_evalue+e_minus, math.e)-math.log(min_outgroup_evalue+e_minus, math.e), '.2f')
        try :
            outg_pct = format(len(set(outgroup_species))/(len(set(outgroup_species))+len(set(ingroup_species))), '.4f')
        except :
            continue

        # print(outgroup_dict)
        # print(ingroup_dict)
        # print(min_outgroup_evalue)
        # print(min_ingroup_evalue)
        print('Alien index:', alienindex)
        print('Outg_pct:', outg_pct)

        if float(alienindex) > 45 and float(outg_pct) > 0.90 and len(set(outgroup_species)) > 30 :
            print('This is a HGT event')
            print('Accession_id: %s' % min_outgroup_key)
            print('Evalue: %s' % evalue[min_outgroup_key])
            taxid= getTaxid(min_outgroup_key)
            Evalue = evalue[min_outgroup_key]
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            # print(lineage)
            print(ranks2lineage)
            # taxid2name = ncbi.get_taxid_translator([ranks2lineage['species'], ranks2lineage['phylum'], ranks2lineage['superkingdom']])
            try :
                taxid2name = ncbi.get_taxid_translator([ranks2lineage['species'], ranks2lineage['phylum'], ranks2lineage['superkingdom']])
                print(taxid2name[ranks2lineage['phylum']])

                superkingdom = taxid2name[ranks2lineage['superkingdom']]
                phylum = taxid2name[ranks2lineage['phylum']]
                hgt = 'Yes'
                taxonomy = superkingdom + '/' + phylum

                item = [gene, alienindex, outg_pct, Evalue, min_outgroup_key, taxonomy]
                if superkingdom == 'Bacteria' :
                    tsv_writer.writerow(item)
            except :
                continue

        else:
            # print('This is not a HGT')
            hgt = 'No'

    outfile.close()


if __name__== "__main__":
    main()
