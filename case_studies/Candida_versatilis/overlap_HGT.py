#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import xlrd
import pandas as pd


def main() :
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

if __name__ == "__main__" :
    main()
