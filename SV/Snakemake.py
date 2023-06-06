## Author: Meenu Bhati

## parsing the required files
configfile: "config.yaml"

# Deal with the lsf logfile
import os
os.system("mkdir loglsf")

# Reduced paths
BAMDIR  = config['resources']['data']
OUTDIR  = config['foldout']['out']
REF     = config['resources']['ref']
ENV     = config ['resources']['env']
sample_file = config['resources']['samples']

sample_file = open(sample_file, 'r')
SAMPLES = [value.strip().split('\t')[0] for value in sample_file.readlines()] #only the first column is taken, in case there are additional columns

rule all:
    input:
        step1_call = expand( OUTDIR + "step1_call/{sample}-smoove.genotyped.vcf.gz",sample=SAMPLES),
        step2_merge = OUTDIR + "step2_merge/smoove_merge.sites.vcf.gz",
        step3_genotype = expand( OUTDIR + "step3_genotype/{sample}-smoove.genotyped.vcf.gz",sample=SAMPLES),
        step4_paste = OUTDIR + "step4_paste/smoove_combined.smoove.square.vcf.gz"

        
rule smoove_call:
    input:
        ref = REF + "ARS-UCD1.2_Btau5.0.1Y.fa",
        bam = BAMDIR + "{sample}.bam",
        bai = BAMDIR + "{sample}.bam.bai"
    output:
        vcf = OUTDIR + "step1_call/{sample}-smoove.genotyped.vcf.gz",
        idx = OUTDIR + "step1_call/{sample}-smoove.genotyped.vcf.gz.csi" 
    params:
        contigs= "~^NKL",
        sample = "{sample}",
        out = OUTDIR + "step1_call"
    conda:
        ENV + "pysmmove.yaml"
    shell:
        """  
        cd $TMPDIR
        
        pwd
        
        cp {input.bam} .
        cp {input.bai} .
        
        smoove call \
        --name {params.sample} \
        --outdir . \
        --fasta {input.ref} \
        --duphold \
        --genotype \
        -p 1 \
        --excludechroms {params.contigs} \
        --exclude /cluster/work/SV/Gap_ARS-UCD_assembly_autosome_50bp_updown.bed \
        {wildcards.sample}.bam
        
        mkdir -p {params.out}
        
        cp $TMPDIR/{wildcards.sample}-smoove.genotyped.vcf.gz* {params.out}
        
        
        """

rule smoove_merge:
    input:
        ref = REF + "ARS-UCD1.2_Btau5.0.1Y.fa",
        vcfs = expand(OUTDIR + "step1_call/{sample}-smoove.genotyped.vcf.gz",sample=SAMPLES),
        idx = expand(OUTDIR + "step1_call/{sample}-smoove.genotyped.vcf.gz.csi",sample=SAMPLES)
    output:
        vcf = OUTDIR + "step2_merge/smoove_merge.sites.vcf.gz",
        counts = OUTDIR + "step2_merge/smoove_merge.smoove-counts.html"
    params:
        out = OUTDIR + "step2_merge",
        name = "smoove_merge"
    conda:
        ENV + "pysmmove.yaml"
    shell:
        """
        smoove merge --name {params.name} \
        --outdir {params.out} \
        -f {input.ref} \
        {input.vcfs}
        """
        
rule smoove_genotype:
    input:
        ref = REF + "ARS-UCD1.2_Btau5.0.1Y.fa",
        vcf = rules.smoove_merge.output.vcf,
        bam = BAMDIR + "{sample}.bam",
        bai = BAMDIR + "{sample}.bam.bai"
    output:
        vcf = OUTDIR + "step3_genotype/{sample}-smoove.genotyped.vcf.gz",
        idx = OUTDIR + "step3_genotype/{sample}-smoove.genotyped.vcf.gz.csi"
    params:
        out = OUTDIR + "step3_genotype",
        sample = "{sample}"
    conda:
        ENV + "pysmmove.yaml"
    shell:
        """
        cd $TMPDIR
        
        pwd
        
        cp {input.vcf} .
        cp {input.bam} .
        cp {input.bai} .
        
        
        smoove genotype  \
        --name {params.sample} \
        --outdir . \
        --fasta {input.ref} \
        --duphold \
        --removepr \
        -p 1 \
        --vcf smoove_merge.sites.vcf.gz \
        {wildcards.sample}.bam
        
        
        mkdir -p {params.out}
        cp $TMPDIR/{wildcards.sample}-smoove.genotyped.vcf.gz* {params.out}
        
        """       

rule smoove_paste:
    input:
        vcf = expand(OUTDIR + "step3_genotype/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES)
    output:
        out = OUTDIR + "step4_paste/smoove_combined.smoove.square.vcf.gz"
    conda:
        ENV + "pysmmove.yaml"
    params:
        out = OUTDIR + "step4_paste",
        name = "smoove_combined"
    shell:
        """
        smoove paste \
        --name {params.name} \
        --outdir {params.out} {input.vcf}
        """
