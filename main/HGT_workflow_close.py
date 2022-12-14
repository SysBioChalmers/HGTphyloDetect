#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import math
import csv
import warnings 
from Bio import BiopythonWarning
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
        accession_bitscore[accession] = float(blast.strip('\n').split("\t")[3])
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

def main(bitscore_parameter=100, HGTIndex=0.5, out_pct=0.8) :
    name = sys.argv[1:][0]
    if len(sys.argv) == 2 :
        pass
    else :
        bitscore_parameter = float(sys.argv[2].split("=")[1])
        HGTIndex = float(sys.argv[3].split("=")[1])
        out_pct = float(sys.argv[4].split("=")[1])
    genes = list()
    geneSeq = dict()
    HGT = list()
    warnings.simplefilter('ignore', BiopythonWarning)

    with open(name, 'r') as handleGene :
        for record in SeqIO.parse(handleGene, "fasta") :
            gene = str(record.id)
            sequence = str(record.seq)
            geneSeq[gene] = sequence
            genes.append(gene)

    n=0
    for gene in genes :
        n += 1
        print("This is gene %d------------------" %(n))
        print(gene)

        if os.path.exists("./blastp_files/%s.txt" % gene) :
            accession_number,accession_bitscore = parse_NCBI(gene)
            print('Yes, blast file already exists, nice!')
        else :
            # Need to install blast!
            with open('./%s.fasta' % gene, 'w') as outfile :
                outfile.write('>'+gene+'\n')
                outfile.write(geneSeq[gene]+'\n')
            if not os.path.exists("./blastp_files/") :
                os.system('mkdir ./blastp_files')
            myCmd = "blastp -db nr -remote -query ./%s.fasta -max_target_seqs 250 -task 'blastp-fast' -outfmt '7 qacc sacc evalue bitscore length pident' -out ./blastp_files/%s.txt" %(gene, gene)
            os.system(myCmd)
            os.remove('./%s.fasta' % gene)
            accession_number,accession_bitscore = parse_NCBI(gene)

        ncbi = NCBITaxa()
        recipient_accession = list()
        outgroup_accession = list()
        recipient_accession_bitscore = dict()
        outgroup_accession_bitscore = dict()
        recipient_species = list()
        outgroup_species = list()

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
                taxid = getTaxid(accession)
                print('Taxid by BLAST:', taxid)
            except :
                continue

            try :
                lineage = ncbi.get_lineage(taxid)
                lineage2ranks = ncbi.get_rank(lineage)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
                # print(ranks2lineage)

                taxid2name = ncbi.get_taxid_translator(lineage)
                taxonomy_alignment = ranks2lineage
            except :
                print('Warning: %s taxid not found!' % str(taxid))
                continue

            try :
                if taxonomy_alignment['subphylum'] == gene_subphylum :
                    recipient_accession.append(accession)
                    recipient_species.append(taxonomy_alignment['species'])

                if taxonomy_alignment['kingdom'] == gene_kingdom and taxonomy_alignment['subphylum'] != gene_subphylum :
                    outgroup_accession.append(accession)
                    outgroup_species.append(taxonomy_alignment['species'])
            except :
                continue

        for accession_id in recipient_accession :
            recipient_accession_bitscore[accession_id] = accession_bitscore[accession_id]

        for accession_id in outgroup_accession :
            outgroup_accession_bitscore[accession_id] = accession_bitscore[accession_id]

        if recipient_accession_bitscore :
            max_recipient_organism_accession_key = max(recipient_accession_bitscore,key=recipient_accession_bitscore.get)
            max_recipient_organism_bitscore = recipient_accession_bitscore[max_recipient_organism_accession_key]

        if outgroup_accession_bitscore :
            max_outgroup_accession_key = max(outgroup_accession_bitscore,key=outgroup_accession_bitscore.get)
            max_outgroup_bitscore = outgroup_accession_bitscore[max_outgroup_accession_key]
            if max_outgroup_accession_key :
                max_taxid = getTaxid(max_outgroup_accession_key)
                max_lineage = ncbi.get_lineage(max_taxid)
                max_lineage2ranks = ncbi.get_rank(max_lineage)
                max_ranks2lineage = dict((rank, taxid) for (taxid, rank) in max_lineage2ranks.items())
                try :
                    max_taxid2name = ncbi.get_taxid_translator([max_ranks2lineage['kingdom'], max_ranks2lineage['phylum'], max_ranks2lineage['subphylum'], max_ranks2lineage['species']])
                except :
                    continue

                # print(gene)
                # print(max_recipient_organism_bitscore)
                # print(max_outgroup_bitscore)
                # print(max_taxid2name)

                if recipient_species :
                    recipient_species_number = len(set(recipient_species))
                if outgroup_species :
                    outgroup_species_number = len(set(outgroup_species))

                HGT_index = format(max_outgroup_bitscore/max_recipient_organism_bitscore, '.4f')
                Outg_pct = format(outgroup_species_number/(outgroup_species_number+recipient_species_number), '.4f')

                print('HGT index: %s' % str(HGT_index))
                print('Out_pct: %s' % str(Outg_pct))
                if max_outgroup_bitscore>=bitscore_parameter and float(HGT_index)>=HGTIndex and float(Outg_pct)>=out_pct :
                    print('This is a HGT event')
                    taxonomy = max_taxid2name[max_ranks2lineage['kingdom']] + '/' + max_taxid2name[max_ranks2lineage['subphylum']]
                    item = [gene, max_outgroup_bitscore, Outg_pct, HGT_index, taxonomy]
                    HGT.append(item)
                else :
                    print('This is not a HGT event')
                    item = [gene, max_outgroup_bitscore, Outg_pct, HGT_index, 'No']
                    HGT.append(item)

    outfile = open("./output_close_HGT.tsv", "wt")
    tsv_writer = csv.writer(outfile, delimiter="\t")
    column = ['Gene/Protein', 'Bitscore', 'Out_pct', 'HGT index', 'Donor taxonomy']
    tsv_writer.writerow(column)
    for HGT_info in HGT :
        tsv_writer.writerow(HGT_info)
    outfile.close()


if __name__== "__main__":
    main()
