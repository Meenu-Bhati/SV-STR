#!/usr/bin/env python
# coding: utf-8

file_name="Bos_taurus.ARS-UCD1.2.104.chr.gtf3"
file=open(file_name,'rt')
contents = file.readlines()
All_gene=dict()

for i in range(len(contents)):
    fline=contents[i].rstrip().split('\t')
    All_gene[i] = [0,0,0,0,0,0,0,0,0,0,0,0]
    All_gene[i][0:6]=fline[0],fline[2],fline[3],fline[4],fline[6],fline[7]
    gene_ids=fline[8].split(";")
    for string in gene_ids:
        ids=string.split()
        if 'gene_id' in ids: 
            All_gene[i][6]=ids[1].replace('"', '')
        elif 'gene_name' in ids:
            All_gene[i][7]=ids[1].replace('"', '')    
        elif 'gene_biotype' in ids:
            All_gene[i][8]=ids[1].replace('"', '')
        elif 'exon_id' in ids:
            All_gene[i][9]=ids[1].replace('"', '')
        elif 'transcript_id' in ids:
            All_gene[i][10]=ids[1].replace('"', '')
        elif 'protein_id'in ids:
            All_gene[i][11]=ids[1].replace('"', '')

            ### write
path="/cluster/work/Annotation/"

All_gene_file = '%s%s%s' % (path, "ARS_UCD104", "_All_104_newcomplete.txt" )
All_gene_out = open (All_gene_file,  "w")

for key, value in All_gene.items():
    All_gene_out.write('%s\t%s\n' % (key,'\t'.join(str(v) for v in value)))

All_gene_out.close ()
