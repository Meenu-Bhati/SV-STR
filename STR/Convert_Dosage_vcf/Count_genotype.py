#!/usr/bin/env python
# coding: utf-8 
## Author: Meenu Bhati  
# In[ ]:
    

#!/usr/bin/env python
# coding: utf-8

# In[4]:

import sys, os, string, time, gzip, re
import numpy as np
import pandas as pd

try:
    path = sys.argv[1]
    filename = sys.argv[2]
    prefix = sys.argv[3]
    indi = sys.argv[4]
except:
    sys.exit()

# In[5]:

vcf_gz='%s%s' % (path,filename)
file=gzip.open(vcf_gz,'rt')
All_pos=dict()
All_GT=dict()

Keep_pos=dict()

for line in file:
    if not line.startswith('##'):
        fline=line.rstrip().split('\t')
        if fline[0] == "#CHROM":
            fids=fline[9:]
        else:
            mypos="\t".join( [fline[0],fline[1]])
            fgts=fline[9:]
            fgts=[el for el in fgts]
            genotype=dict(zip(fids,fgts))
            missing = str(fgts.count("."))
            total = str(len(fgts))
            genotype_count= "\t".join( [total,missing])
            counts = {item:fgts.count(item) for item in fgts}
            All_pos[mypos] = genotype_count
            All_GT[mypos] = counts
            #unwanted_set={'.', '0'}
            unwanted_set={'.'}
            actual_set=set(counts)
            ans=set(actual_set)-unwanted_set
            miss=fgts.count(".")
            if len(ans) > 1 and miss < indi :
                Keep_pos[mypos] = genotype_count

## keep pos is to keep those position where DS is atleast 2 genotypes other than .
                
pos_out_file = '%s%s%s' % (path, prefix, "_genotype_count.txt"  )
pos_out = open (pos_out_file,  "w")

for pos in All_pos:
    pos_out.write (f"{pos}\t{All_pos.get(pos)} \n")
pos_out.close ()

pos_out_all = '%s%s%s' % (path, prefix, "_All_geno_count.txt"  )
pos_out2 = open (pos_out_all,  "w")

for val in All_GT:
    pos_out2.write (f"{val}\t{All_GT.get(val)} \n")
pos_out2.close ()

### keep the sites with DS >=1 other than 0 and .  and missing sites less than 15 

pos_out_keep = '%s%s%s' % (path, prefix, "_pos_keep.txt"  )
pos_out3 = open (pos_out_keep,  "w")

for val in Keep_pos:
    pos_out3.write (f"{val}\t{Keep_pos.get(val)} \n")
pos_out2.close ()



