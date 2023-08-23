#!/usr/bin/env python
# coding: utf-8
## Author: Meenu Bhati
# In[2]:


import argparse
import gzip
import math
import numpy as np
import os
from collections import Counter
import pandas as pd
import random
import shutil
import sys
import tabix
import pysam
from itertools import chain


def GetKeepGTs(DS, MINGT):
    all_ds = Counter(chain(DS))
    keepgt = [item for item in all_ds if None not in item and all_ds [item] >= MINGT]
    return keepgt

parser = argparse.ArgumentParser( description='Change DS to none if no of DS is less than MIN GT')
parser.add_argument("-v", "--invcf", help="SV/SNP input VCF_file", nargs=1, type=str)
parser.add_argument("-o", "--outvcf", help="SV/SNP output VCF_file", nargs=1,type=str)
parser.add_argument("--mingt", help="Remove genotypes with fewer than this many samples", type=int, default=1)

if __name__ == '__main__':
    args = parser.parse_args()
    input_vcf = args.invcf[0]
    out_vcf = args.outvcf[0]
    MINGT = args.mingt

myvcf = pysam.VariantFile(input_vcf, "r")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile( out_vcf, "w", header=myvcf.header)

vcf_out = pysam.VariantFile( out_vcf, "w", header=myvcf.header)

myvcf = pysam.VariantFile(input_vcf, "r")

with open(out_vcf, "a") as out:
    for variant in myvcf:
        DS = [value['DS'] for value in variant.samples.values()]
        keepgt = GetKeepGTs (DS, MINGT)
        for sample in variant.samples:
            if variant.samples[sample]['DS'] not in keepgt:
                variant.samples[sample]['DS'] = None
        out.write(str(variant))

