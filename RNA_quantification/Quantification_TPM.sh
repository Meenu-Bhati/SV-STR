### check mutiqc on final bam

### 1. after alignment

final_bam="dedup_alignment"
python -m multiqc ${final_bam}/* --interactive

### 2. Gene level quantification

##path to dir
input="path/dedup_alignment"

samples="breed_coverage.txt" ### here first coloum is sample id of bam file e.g. sample1.bam

quant_path="path/quantification"

mkdir -p ${output_path}/logs


awk '{ print $1}' $samples | while read sample
do
QTLtools quan --bam ${input}/${sample}/${sample}.bam --gtf ensembl_104/Bos_taurus.ARS-UCD1.2.104.chr.gtf3.gz --sample ${sample} --out-prefix ${quant_path}/${sample} --filter-mapping-quality 60 --check-proper-pairing --filter-failed-qc --tpm
done


### 3. combine in one file

# create one file that includes all samples

quant_path="/cluster/work/pausch/meenu/SV/expression/cis_eQTL_new/quantification"


cat ${quant_path}/sample1.*.gene.tpm.bed  | cut -f1-6 > ${quant_path}/anno_info

awk '{ print $1}' $samples | while read sample
do
cat ${quant_path}/${sample}.*.gene.tpm.bed | cut -f7 > ${quant_path}/${sample}_temp_quant 
done

paste -d "\t" ${quant_path}/anno_info ${quant_path}/*_temp_quant > ${quant_path}/FINAL_gene_tpm.tsv

rm ${quant_path}/*_temp_quant
sed -i  's/#//g' FINAL_gene_tpm.tsv


### 4. now normalize data

Rscript TPM_normalization.r


### 5. bgzip the data

bgzip FINAL_gene_TPM_filter_quan_inv_allgenes.bed


