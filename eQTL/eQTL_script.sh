### 1.  nominal eQTL 


## path to files
module load gcc/8.2.0 gsl

expr="FINAL_gene_TPM_filter_quan_inv_allgenes.bed.gz"
cov="All_final_PEER_covariates.txt"

qtltools=qtltools/bin/QTLtools


out="nominal"

## this is combined SV and STR vcf 
vcf_file=STRmac1_filter_SV_maf0.05_geno80_75ids.norm_sorted.vcf.gz

## Run command 
${qtltools}  cis --vcf ${vcf_file} --bed ${expr} \
--std-err  --normal --nominal 0.01  --cov ${cov} --out  ${out}/joint_sv_str_nominal --window 1000000



### 2. Permutation eQTL


qtltools=/cluster/work/pausch/alex/software/qtltools/bin/QTLtools

out_perm="permutation" ## dir to save results 

${qtltools} cis --vcf ${vcf_file} --bed ${expr} --cov ${cov} \
--std-err --normal --permute 1000 --out ${out_perm}/joint_str_sv_perm --window 1000000 

## header
${qtltools} cis --vcf ${vcf_file} --bed ${expr} --cov ${cov} \
--std-err --normal --permute 1000 --out ${out_perm}/header_perm --window 1000000 --chunk 0 1 


## gzip
module load htslib

cat ${out_perm}/header_perm ${out_perm}/joint_str_sv_perm > ${out_perm}/joint_str_sv_perm.txt

gzip  ${out_perm}/joint_str_sv_perm.txt

## fdr 
module load gcc/4.8.5 r/4.1.3

#qtltools_runFDR_cis.R script is provided by QTLtools
Rscript qtltools_runFDR_cis.R \
${out_perm}/joint_str_sv_perm.txt.gz 0.05 \
${out_perm}/joint_str_sv_perm_fdr > log_rscript_perm



### 3. Conditional eQTL


out_condi="conditional"

${qtltools} cis --vcf ${vcf_file} \
--bed ${expr} --cov ${cov} --mapping ${out_perm}/joint_str_sv_perm_fdr.thresholds.txt \
--std-err --out ${out_condi}/joint_str_sv_condi --window 1000000  

### header

${qtltools} cis --vcf ${vcf_file} \
--bed ${expr} --cov ${cov} --mapping ${out_perm}/joint_str_sv_perm_fdr.thresholds.txt  \
--std-err  --out ${out_condi}/header_condi --window 1000000 --chunk 0 1

cat ${out_condi}/header_condi ${out_condi}/joint_str_sv_condi > ${out_condi}/joint_str_sv_condi.txt

