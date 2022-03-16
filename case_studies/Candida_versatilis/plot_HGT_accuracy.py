#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import math
import xlrd
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def get_genes() :
    # get the protein sequence accoding to protein sequence id
    with open("./Candida_versatilis.fasta", "r") as handleGene :
        genes = list()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            genes.append(str(record.id))
        print("The gene number of this organism is:", len(genes))
    return genes

def preprocess() :
    genes_all = get_genes()
    print(len(genes_all))  # 5475
    worksheet = xlrd.open_workbook(u"./HGT_case.xlsx")
    sheet_names = worksheet.sheet_names()
    # print(sheet_names)
    sheet = worksheet.sheet_by_name(sheet_names[0])
    rows = sheet.nrows
    # print(rows)
    genes = list()
    for i in range(1,rows) :
        cell = sheet.cell_value(i,1)
        # print(cell)
        if cell.startswith('Candida_versatilis') :
            genes.append(cell)
    HGT_genes_cell = genes
    print(len(HGT_genes_cell))  # 169
    print(HGT_genes_cell[-10:])

    with open("./candida_versatilis_HGT.tsv", "r") as file :
        lines = file.readlines()[1:]

    HGT_genes_us = [line.strip().split('\t')[0] for line in lines]
    print(len(HGT_genes_us))  # 228
    print(HGT_genes_us[:10])

    HGT_interac = [HGT for HGT in HGT_genes_us if HGT in HGT_genes_cell]

    print(len(HGT_interac))  # 148  
    HGT_detection_ratio = '%.4f' % (len(HGT_interac)/len(HGT_genes_cell))
    print('The percentage of HGT detection number in all HGT events:', HGT_detection_ratio)  # 87.57%

    return genes_all, HGT_genes_cell, HGT_genes_us

def calculation() :
    genes_all, HGT_genes_cell, HGT_genes_us = preprocess()
    genes_HGT = dict()
    for gene in genes_all :
        if gene in HGT_genes_cell :
            genes_HGT[gene] = 'HGT'
        else :
            genes_HGT[gene] = 'Non_HGT'

    tp, tn, fp, fn = 0, 0, 0, 0
    for gene in genes_all :
        if genes_HGT[gene] == 'HGT' :
            if gene in HGT_genes_us :
                tp += 1
            else :
                fn += 1
        else :
            if gene not in HGT_genes_us :
                tn += 1
            else :
                fp += 1
    print(tp, tn, fp, fn)  # 148 5226 80 21
    return tp, tn, fp, fn

def main() :
    tp, tn, fp, fn = calculation()
    my_metrics = {
        'Sensitivity': 'NA',
        'Specificity': 'NA',
        'Accuracy': 'NA',
        'MCC': 'NA',
        'Recall': 'NA',
        'Precision': 'NA',
        'F1-score': 'NA',
    }

    my_metrics['Sensitivity'] = '%.4f' % (tp / (tp + fn)) if (tp + fn) != 0 else 'NA'
    my_metrics['Specificity'] = '%.4f' % (tn / (fp + tn)) if (fp + tn) != 0 else 'NA'
    my_metrics['Accuracy'] = '%.4f' % ((tp + tn) / (tp + fn + tn + fp))
    my_metrics['MCC'] = '%.4f' % ((tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))) if (tp + fp) * (
        tp + fn) * (tn + fp) * (tn + fn) != 0 else 'NA'
    my_metrics['Precision'] = '%.4f' % (tp / (tp + fp)) if (tp + fp) != 0 else 'NA'
    my_metrics['Recall'] = my_metrics['Sensitivity']
    my_metrics['F1-score'] = '%.4f' % (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) != 0 else 'NA'

    print(my_metrics)
    # {'Sensitivity': '0.8757', 'Specificity': '0.9849', 'Accuracy': '0.9816', 'MCC': '0.7451', 'Recall': '0.8757', 'Precision': '0.6491', 'F1-score': '0.7456'}
    # return my_metrics

    metrics = ['Accuracy', 'Sensitivity', 'Specificity', 'MCC', 'F1-score'] 
    percentage = [float(my_metrics['Accuracy'])*100, float(my_metrics['Sensitivity'])*100, float(my_metrics['Specificity'])*100, float(my_metrics['MCC'])*100, float(my_metrics['F1-score'])*100]
    # print(type(my_metrics['Accuracy']))  # str

    plt.figure(figsize=(4.0,3.0))
    plt.bar(range(len(percentage)), percentage, tick_label=metrics, width=0.5, alpha=0.8, color='pink', edgecolor='r')

    # plt.ylim([0,1.0])
    # plt.yticks([0,0.20,0.40,0.60,0.80,1.00])

    plt.ylim([0,100])
    plt.yticks([0,20,40,60,80,100])

    # plt.xlabel('Metrics',fontsize=12)
    plt.ylabel("Score (%)", fontsize=12)
    plt.xticks(fontsize=12, rotation=30, ha='right')
    plt.yticks(fontsize=12)

    plt.savefig("./metrics_score.pdf", dpi=400, bbox_inches='tight')


if __name__ == "__main__" :
    main()
