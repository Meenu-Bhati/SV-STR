## Author: Meenu Bhati

## parsing the required files
configfile: "config.yaml"

# Deal with the lsf logfile
import os
os.system("mkdir log_folder")

# Reduced paths
reference =  config['path']['ref']
ref_STR = config ['path']['ref_str']
output_path   =  config['path']['output']
bamfiles = config ['path']['bamfiles']
filter_out = config['path']['outfilter']

#tools
HIPSTR = config['tools']['hipstr']
ENV     = config ['tools']['env']
SCRIPT = config ['tools']['script']
VEP   = config ['tools']['vepimage']
BCFTOOLS = config ['tools']['bcftools']
## chr
chrs = list(range(1,30)) + ['X'] ## same as in ref dict

rule all:
    input:
        hip_out= expand ( output_path + "{chr}_str_gts.vcf.gz", chr=chrs),
        filt_out= expand ( filter_out + "/{chr}_str_gts_filter.vcf.gz",chr=chrs),
        concatvcf= filter_out + "/Autosomes_str_filter.vcf.gz"
        #anno= filter_out + "/Autosomes_str_filter.vcf.vep.gz"
        
rule hipstr_run:
    input: 
        ref_str = ref_STR   + "ARS_UCD.hipstr_reference_All_info.bed",
        bamfiles =  bamfiles + "bamfiles_183_noBGI.txt"   
    output:
        vcf = output_path + "{chr}_str_gts.vcf.gz"
    params:
        Chr = "{chr}"     
    shell:
        HIPSTR + " --bam-files {input.bamfiles}" + " --fasta " + reference + " --regions {input.ref_str} --chrom {params.Chr} --def-stutter-model --min-reads 10 --log log_folder/{params.Chr}_log.txt --str-vcf {output.vcf}"
             
rule filter_call:
    input:
        vcf = output_path + "{chr}_str_gts.vcf.gz"
    output:
        vcf = filter_out +  "/{chr}_str_gts_filter.vcf.gz"
    conda: 
        ENV + "/py2.yaml"
    shell:
        """
        python Filter_STR.py  --vcf {input.vcf} \
        --min-call-qual 0.8 \
        --max-call-flank-indel  0.20 \
        --min-loc-depth 5  \
        > {output.vcf} 
                
        """
        
rule concat_vcf:
    input:
        calls = expand ( filter_out +  "/{chr}_str_gts_filter.vcf.gz", chr=chrs)
    output:
        multiext(filter_out +"/Autosomes_str_filter.vcf",".gz",".gz.tbi")
    shell:
        """
        module load htslib/1.12
        
        bcftools concat -o {output[0]} -Oz {input.calls}
        
        tabix -p vcf {output[0]}
        
        """        

