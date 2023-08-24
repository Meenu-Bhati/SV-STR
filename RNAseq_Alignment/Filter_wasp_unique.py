#!/usr/bin/env python
# coding: utf-8
## Author: Meenu Bhati
# In[ ]:


## script to filter out reads which did not pass WASP filtering

##vW tag to the SAM output: vW:i:1 means alignment passed WASP  filltering, and all other values mean it did not pass:

import sys, os, string, time, gzip, re
import pysam

try:
    infile = sys.argv[1]
    outfile = sys.argv[2]
except:
    print("provide infile bam path and outfile bam path with filename ")
    sys.exit()
##  out="new_vW.bam" example 

samfile=pysam.AlignmentFile(infile, "rb")
    
with pysam.AlignmentFile(outfile, "wb", template=samfile) as outf:
    for x in samfile:
        if x.mapping_quality == 60:
            tags=x.tags
            output = [i[1] for i in tags if i[0] == 'vW']
            if len(output) == 0:
                outf.write(x)
            elif 1 in output:
                outf.write(x)

print("Filtering done")

