## Annotate cluster ids

## 1. first prepare file of cluster ids with chr, start, end, clusterid, strand

awk '{ print $1}' Top_conditional_joint_all_condi | sort | uniq | awk -v OFS="\t" 'NR>1{split($1,arr,":" ); split(arr[4], sign, "_"); print arr[1],arr[2],arr[3],$1,sign[3] }' > Top_conditional_joint_all_condi_cluster_ids
sed -i '1s/^/chr\tstart\tend\tpid\tstrand\n/' Top_conditional_joint_all_condi_cluster_ids

## 2. Annotate

Rscript Annontate_introncluster.r \
--exon Complete_exon_format.txt --intron Top_conditional_str_condi_cluster_ids --topfile Top_conditional_str_condi_header.txt --outfile Top_conditional_str_condi_all_info.txt
