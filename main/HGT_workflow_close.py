#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import math
import csv
from Bio import SeqIO
from Bio import Entrez
from ete3 import NCBITaxa


def parse_NCBI(gene):
    with open("./blastp_files/%s.txt" % gene, "r") as filename :
        data = filename.readlines()
        for line in data :
            if line.strip("\n").endswith("found") :
                index = data.index(line)

        blast_results = data[index+1:-1]

    accession_number = list()
    accession_bitscore = dict()
    for blast in blast_results :
        accession = blast.strip("\n").split("\t")[1]
        accession_number.append(accession)
        accession_bitscore[accession] = float(blast.strip('\n').split("\t")[-2])
    return accession_number, accession_bitscore

def getTaxid(accession):
    # Retrieving data in the GenBank using only the GenBank code accession in biopython
    # https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    # https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    # https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    Entrez.email = "abcd@ncbi.org"

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

    return taxid

def main() :
    name = sys.argv[1:][0]
    genes = list()
    HGT = list()

    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            genes.append(gene)

    n=0
    for gene in genes :
        n += 1
        print("This is gene %d------------------" %(n))
        print(gene)

        if os.path.exists("./blastp_files/%s.txt" % gene) :
            accession_number,evalue = parse_NCBI(gene)
            print('Yes, blast file already exists, nice!')
        else :
            # Need to install blast!
            myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue length pident' -out ./blastp_files/%s.txt" %('input', gene)
            os.system(myCmd)

        ncbi = NCBITaxa()
        recipient_organism = list()
        outgroup_accession = list()
        recipient_accession_bitscore = dict()
        outgroup_accession_bioscore = dict()
        recipient_species = list()
        outgroup_species = list()

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
                    recipient_organism.append(accession)
                    recipient_species.append(taxonomy_alignment['species'])

                if taxonomy_alignment['kingdom'] == 4751 and taxonomy_alignment['subphylum'] != 147537 :
                    outgroup_accession.append(accession)
                    outgroup_species.append(taxonomy_alignment['species'])
            except :
                continue

        for accession_id in recipient_organism :
            recipient_accession_bitscore[accession_id] = accession_bitscore[accession_id]

        for accession_id in outgroup_accession :
            outgroup_accession_bioscore[accession_id] = accession_bitscore[accession_id]

        if recipient_accession_bitscore :
            max_recipient_organism_accession_key = max(recipient_accession_bitscore,key=recipient_accession_bitscore.get)
            max_recipient_organism_bitscore = recipient_accession_bitscore[max_recipient_organism_accession_key]

        if outgroup_accession_bioscore :
            max_outgroup_accession_key = max(outgroup_accession_bioscore,key=outgroup_accession_bioscore.get)
            max_other_species_bitscore = outgroup_accession_bioscore[max_outgroup_accession_key]
            if max_outgroup_accession_key :
                max_taxid = getTaxid(max_outgroup_accession_key)
                max_lineage = ncbi.get_lineage(max_taxid)
                max_lineage2ranks = ncbi.get_rank(max_lineage)
                max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
                try :
                    max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage['kingdom'], max_ranks2lineage['phylum'], max_ranks2lineage['subphylum'], max_ranks2lineage['species']])
                except :
                    continue

                print(gene)
                print(max_recipient_organism_bitscore)
                print(max_other_species_bitscore)
                print(max_taxid2name)

                if recipient_species :
                    recipient_species_number = len(set(recipient_species))
                if outgroup_species :
                    outgroup_species_number = len(set(outgroup_species))

                HGT_index = format(max_other_species_bitscore/max_recipient_organism_bitscore, '.4f')
                Outg_pct = format(outgroup_species_number/(outgroup_species_number+recipient_species_number), '.4f')

                print(HGT_index)
                print(recipient_species_number)
                print(outgroup_species_number)
                print(Outg_pct) 
                if max_other_species_bitscore>100 and float(HGT_index)>=0.5 and float(Outg_pct)>=0.8 :
                    taxonomy = max_taxid2name[max_ranks2lineage['kingdom']] + '/' + max_taxid2name[max_ranks2lineage['subphylum']]
                    item = [gene, max_other_species_bitscore, Outg_pct, HGT_index, taxonomy]
                    HGT.append(item)

    outfile = open("./output.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    column = ['Gene', 'Bitscore', 'Outg_pct', 'HGT index', 'Donor taxonomy']
    tsv_writer.writerow(column)
    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)
    outfile.close()


if __name__== "__main__":
    main()
