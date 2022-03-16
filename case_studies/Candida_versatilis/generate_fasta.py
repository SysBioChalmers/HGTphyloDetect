#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# This python script is to generate the fasta file for Candida_versatilis

import json
from Bio import SeqIO


def get_refSeq(organism) :
    # get the protein sequence accoding to protein sequence id
    with open("/Users/leyu/Documents/Le/Data/orthomcl_output/343taxa_proteins.fasta", "r") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            if record.id.startswith(organism) :
                proteinSeq[record.id] = str(record.seq)
        print("The protein number of this organism is:", len(proteinSeq))
        # for key in proteinSeq.keys() :
        #     print(key)
        #     print(proteinSeq[key])
    return proteinSeq

def main() :
    proteinSeq = get_refSeq("Candida_versatilis")
    with open("./Candida_versatilis.fasta", "w") as f :
        for proteinId in proteinSeq.keys() :
            f.write(">%s\n" % (proteinId))
            f.write("%s\n" % (proteinSeq[proteinId]))


if __name__ == "__main__" :
    main()
