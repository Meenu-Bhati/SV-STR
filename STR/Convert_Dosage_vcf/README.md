## Steps for converting VCF file into dosage VCF format, compatible with QTLtools.

### 1. Add dosage score (DS) in the format field of vcf
Dosage score are sum of GB tag

```
module load bcftools
module load htslib
        
python Add_dosage_vcf.py  --invcf  ${dir}/str.vcf.gz --outvcf str_DSinfo.vcf
bgzip  str_DSinfo.vcf
tabix -p vcf str_DSinfo.vcf.gz
```
### 2. Keep only DS field in the vcf

```
bcftools annotate -x ^FORMAT/DS -O z -o str_DS.vcf.gz str_DSinfo.vcf.gz
```

### 3. Replace Alt allele with "STR" and change DS type to float in vcf header 
It is necessary; otherwise, multiallelic sites will be removed from QTLtools. 

```
module load htslib
        
zcat str_DS.vcf.gz |  awk '$1 ~ /^#/ {OFS="\t"; print $0;next} {$5="STR"; print $0}' | sed '2244s/.*/##FORMAT=<ID=DS,Number=A,Type=Float,Description=estimated ALT dose/' > str_DS_alt.vcf


## 2244 is line number of ##FORMAT field of DS in my vcf ## Check your file and replace this number with correct

```

### 4. Convert DS to missing(.) if it is observed in less than n number of animals 
This step aims to minimize the impact of outliers

```
## in this study n=2
python Remove_genotype_with_less_number.py -v str_DS_alt.vcf -o str_DS_alt_minGT.vcf --mingt n
bgzip  str_DS_alt_minGT.vcf
tabix -p str_DS_alt_minGT.vcf.gz

```
### 5. Filter STR loci  atleast 80% genotyped and atleast more than 2 DS scores per loci 

```
## in this study individuals=15
## prefix_output=filter
python Count_genotype.py <path_str_file>  str_DS_alt_minGT.vcf.gz  prefix_output individuals

#The output of pos keep file
#chr pos   N    missing_n
#1	215422	75	10
#1	220167	75	7
#1	220373	75	0

## Keep only sites high sites 
awk -v OFS="\t" '{ print $1,$2}' filter_pos_keep.txt > postion.txt

bcftools view --regions-file postion.txt -O z -o STR_final.vcf.gz  str_DS_alt_minGT.vcf.gz

tabix -p STR_final.vcf.gz

```


