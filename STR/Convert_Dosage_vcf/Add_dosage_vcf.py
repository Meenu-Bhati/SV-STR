#!/usr/bin/env python
# coding: utf-8
# Author:Meenu Bhati
# In[32]:


import sys
import gzip
import vcf
import argparse
import os
import numpy as np
import pysam


parser = argparse.ArgumentParser( description='Add dosage i.e sum of GB in STR vcf file ')
parser.add_argument("-v", "--invcf", help="STR input VCF_file", nargs=1, type=str)
parser.add_argument("-o", "--outvcf", help="STR output VCF_file", nargs=1,type=str)

if __name__ == '__main__':
    args = parser.parse_args()
    input_vcf = args.invcf[0]
    out_vcf = args.outvcf[0]


myvcf = pysam.VariantFile(input_vcf, "r")

myvcf.header.formats.add("DS", ".", "String", "Sum value of GB")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile( out_vcf, "w", header=myvcf.header)

vcf_out = pysam.VariantFile( out_vcf, "w", header=myvcf.header)


with open(out_vcf, "a") as out:
    for variant in myvcf:
        for sample in variant.samples:
            DS_value = '.'
            if "|" in variant.samples[sample]['GB']:
                DS_value = sum(float(i) if len(variant.samples[sample]['GB']) > 1 else np.nan for i in variant.samples[sample]['GB'].split('|'))
                variant.samples[sample]['DS'] = str(DS_value)
            elif "." in variant.samples[sample]['GB']:
                variant.samples[sample]['DS'] = str(".")
        out.write(str(variant))


# In[ ]:




