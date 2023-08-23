#STEP 1: Generate exon-exon junction counts with regtools

## path to dir
regtools=regtools/build/regtools
samples="sample_id.txt" ### sampleid or prefix of each bam file in first coloumn.. eg. sample1 in sample1.bam 

inputpath="Alignement/wasp_filter" ## wasp filtered bam files
outputpath="Reg_output"

mkdir -p ${outputpath}
awk '{ print $1}' $samples | while read sample
do
${regtools} junctions extract -a 8 -m 50 -M 500000 -s 1 ${inputpath}/${sample}/${sample}_filter.bam -o ${outputpath}/${sample}.regtools.junc.txt
done



## Step 2: Merge all regtools outout 

### bgzip above files

awk '{ print $1}' sample_id.txt | while read sample; do bgzip ${outputpath}/${sample}.regtools.junc.txt; done

### Create a list containing all paths and file names
awk '{ print $1}' sample_id.txt | while read sample; do echo ${outputpath}/${sample}.regtools.junc.txt; done > All_junc_list.txt



## Step 3 Generate intron excision ratios with LeafCutter

module load r/4.0.2 
module load gsl
module load gcc/8.2.0

## things to do 
# create collapsed transcripts per gene
https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model

## create sample map file
awk '{ print $1}' sample_id.txt | while read sample; do \
echo -e ${sample}"\t"${sample}; done > sample_lookup.txt

sed -i '1s/^/sample_id\tparticipant_id\n/' sample_lookup.txt


## run leafcutter
# cluster_prepare_fastqtl.py is modified from https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl/leafcutter/src 
cluster_prepare_fastqtl.py \
/cluster/work/pausch/meenu/SV/expression/cis-sQTL/reg_leaf_bam/All_junc_list.txt \
transcritome_input/Bos_taurus.ARS-UCD1.2.104.genes.exons.txt.gz \
transcritome_input/Bos_taurus.ARS-UCD1.2.104.genes.gtf New \
sample_lookup.txt -o leaf_out/ \
--leafcutter_dir tools/leafcutter/ --num_pcs 10

## this will generate a 29 leaf_out.perind.counts.filtered.gz.qqnorm bed files,  further filter based on exon boundries leaf_out.leafcutter.bed.gz
## Due to the incomplete bovine annotation, we used per chr perind.counts.filtered.gz.qqnorm bed files


## Step 4 Modify output from above step with QTLtools format 

touch final_qqnorm_leafcutter.bed
cp leaf_out_perind.counts.filtered.gz.qqnorm_1 total_qqnorm_leafcutter.bed 

for chr in {2..29}
do
## removing header from each file
tail -n+2 leaf_out_perind.counts.filtered.gz.qqnorm_${chr} >> total_qqnorm_leafcutter.bed 
done

bgzip total_qqnorm_leafcutter.bed 


## now seperate sample coloums and add gene/cluster info in reuired format in the final files
zcat total_qqnorm_leafcutter.bed.gz | awk -v n=5 'OFS="\t" { split($4,arr,":"); for (i=n; i<=NF; i++) printf "%s%s",  $i, (i<NF ? OFS : ORS)}'  > Sample_only_leafcutter_out

zcat total_qqnorm_leafcutter.bed.gz | cut -f1-4 | awk -F "\t" -v OFS="\t" 'NR>1{split($4, arr,":" ); split(arr[4], sign, "_"); print $1,$2,$3,$4, arr[5], sign[3]}'  >  Gene_info_leafcutter_out
## add header in gene inof or 
sed -i '1s/^/#chr\tstart\tend\tpid\tgid\tstrand\n/' Gene_info_leafcutter_out

paste -d "\t" Gene_info_leafcutter_out Sample_only_leafcutter_out > final_leafcutter_qtltools_input.bed

bgzip final_leafcutter_qtltools_input.bed
tabix final_leafcutter_qtltools_input.bed.gz 


