#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Step 1 : BLAST hits were parsed to retrieve associated taxonomic information
# 		using the NCBI's taxonomy database.

# Step 2: Calculate Alien Index (AI) values and out_pct based on the above information.

# Step 3: Output the horizontal gene transfer (HGT) high-throughput identification results.

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
    evalue = dict()
    for blast in blast_results :
        accession = blast.strip("\n").split("\t")[1]
        accession_number.append(accession)
        evalue[accession] = float(blast.strip('\n').split("\t")[2])
    return accession_number, evalue

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
        outgroup = list()
        ingroup = list()
        outgroup_species = list()
        ingroup_species = list()
        outgroup_dict = dict()
        ingroup_dict = dict()
        e_minus = 1e-200

        try :
            gene_taxid = getTaxid(gene)
            gene_lineage = ncbi.get_lineage(gene_taxid)
            gene_lineage2ranks = ncbi.get_rank(gene_lineage)
            gene_ranks2lineage = dict((rank, taxid) for (taxid, rank) in gene_lineage2ranks.items())
            gene_taxonomy_alignment = gene_ranks2lineage
            gene_kingdom = gene_taxonomy_alignment['kingdom']
            gene_subphylum = gene_taxonomy_alignment['subphylum']
            # print(gene_kingdom)
            # print(type(gene_kingdom)) # <class 'int'>
            # print(gene_subphylum)
        except :
            print('Attention: please check the gene accession id!')
            break

        for accession in accession_number[:200] :
            try :
                # print(accession)
                taxid = getTaxid(accession)
                print('Taxid by BLAST:', taxid)
            except :
                continue
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())

            taxonomy_alignment = ranks2lineage
            LINNAEUS_FILTER = ["subphylum","kingdom","superkingdom"]

            try :
                if taxonomy_alignment['kingdom'] != gene_kingdom : 
                    outgroup.append(accession) 
                    outgroup_species.append(taxonomy_alignment['species'])

                if taxonomy_alignment['kingdom'] == gene_kingdom and taxonomy_alignment['subphylum'] != gene_subphylum :
                    ingroup.append(accession)
                    ingroup_species.append(taxonomy_alignment['species'])

            except :
                outgroup.append(accession)
                try :
                    outgroup_species.append(taxonomy_alignment['species'])
                except :
                    continue

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

        print('Alien index: %s' % str(alienindex))

        if float(alienindex) > 45 and float(outg_pct) > 0.90 :
            print('This is a HGT event')
            print('Accession_id: %s' % min_outgroup_key)
            print('Evalue: %s' % evalue[min_outgroup_key])
            taxid= getTaxid(min_outgroup_key)
            Evalue = evalue[min_outgroup_key]
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            # print(lineage)
            # print(ranks2lineage)
            taxid2name = ncbi.get_taxid_translator([ranks2lineage['phylum'], ranks2lineage['class'], ranks2lineage['superkingdom']])
            print(taxid2name[ranks2lineage['phylum']])

            superkingdom = taxid2name[ranks2lineage['superkingdom']]
            phylum = taxid2name[ranks2lineage['phylum']]
            hgt = 'Yes'
            taxonomy = superkingdom + '/' + phylum

            item = [gene, alienindex, Evalue, min_outgroup_key, taxonomy]
            HGT.append(item)

        else:
            # print('This is not a HGT event')
            hgt = 'No'
            HGT.append([gene,alienindex,'No', 'No', 'No'])

    outfile = open("./output.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    column = ['Gene/Protein', 'Alien index', 'E value', 'Donor id', 'Donor taxonomy']
    tsv_writer.writerow(column)
    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)
    outfile.close()


if __name__== "__main__":
    main()
