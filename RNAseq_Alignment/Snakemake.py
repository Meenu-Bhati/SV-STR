#!/usr/bin/env python
# coding: utf-8
## Author: Meenu Bhati
# In[ ]:


# Deal with the lsf logfile
import os
os.system("mkdir log_folder")

configfile: "config.yaml"
sample_file = config['resources']['samples']
sample_file = open(sample_file, 'r')

samples = [value.strip().split(' ')[0] for value in sample_file.readlines()] #only the first column is taken which has bam id
chroms = list(range(1,30)) + ['X','MT']

# Reduced paths
raw_data_path = config['resources']['raw']
assembly = config['resources']['assembly']
# info_file = config["resources"]["info_file"]
fastp_output = config['fold_out']['fastp']
split_fastq = config['fold_out']['split_fastq']
vcf_path = config['resources']['vcf']
alignment = config['fold_out']['Alignment']
sorted_alignment = config['fold_out']['sorted_alignment']
dedup_alignment = config['fold_out']['dedup_alignment']
metrics = config['fold_out']['metrics']
wasp_out = config['fold_out']['wasp_out']
# Tools
FASTP = config['tools']['fastp']
multiqc = config['tools']['multiqc']
FASTQ_SPLITTER = config['tools']['fastq_splitter']
SAMTOOLS = config['tools']['samtools']
SAMBAMBA = config['tools']['sambamba']
LOAD_PICARD_TOOLS = config['tools']['load_picard_tools']
LOAD_JAVA = config['tools']['load_java_jdk']
LOAD_R = config['tools']['load_r']
PICARD_TOOLS = config['tools']['picard_tools']

rule all:
    input:
        multiqc = fastp_output + "multiqc/multiqc_report.html",
        split_fastqs = expand(split_fastq + "{sample}", sample = samples),
        stats_1=expand(sorted_alignment + "{sample}/{sample}.stats",sample=samples),
        #stats_2=expand(dedup_alignment + "{sample}/{sample}.stats",sample=samples),
        index= expand(dedup_alignment + "{sample}/{sample}.bam.bai",sample=samples),
        metrics = expand(metrics + "{sample}/{sample}.alignment_summary_metrics",sample=samples),
        stats = expand(metrics + "{sample}/{sample}_idxstats",sample=samples),
        wasp = expand(wasp_out + "{sample}/{sample}_filter.bam.bai",sample=samples)
        
rule fastp:
    input:
        R1 = raw_data_path + "{sample}_R1.fastq.gz",
        R2 = raw_data_path + "{sample}_R2.fastq.gz"
    output:
        R1_O = temp(fastp_output + "{sample}/{sample}_R1.fastq.gz"),
        R2_O = temp(fastp_output + "{sample}/{sample}_R2.fastq.gz"),
        R_HTML = fastp_output + "{sample}/{sample}_fastp.html",
        R_JSON = fastp_output + "{sample}/{sample}_fastp.json"
    params:
        " --trim_poly_g --trim_poly_x "
    shell:
        FASTP + " -i {input.R1} -o {output.R1_O} -I {input.R2} -O {output.R2_O} -h {output.R_HTML} -j {output.R_JSON} {params} >/dev/null"

rule multiqc:
    input:
        fastp_json = expand(fastp_output + "{sample}/{sample}_fastp.json",sample=samples)
    output:
        fastp_output + "multiqc/multiqc_report.html"
    params:
        json = "-k json"
    shell:
        multiqc + " {params.json} " + fastp_output + "* -o" + fastp_output + "multiqc"

checkpoint split_fastq:
    input:
        R1 = fastp_output + "{sample}/{sample}_R1.fastq.gz",
        R2 = fastp_output + "{sample}/{sample}_R2.fastq.gz"
    output:
        split_fastqs = directory(split_fastq + "{sample}")
    params:
        prefix = split_fastq + "{sample}/{sample}_",
        folder = split_fastq + "{sample}"
    shell:
        "mkdir {params.folder} \n" +
        FASTQ_SPLITTER + " -o {params.prefix} {input.R1} {input.R2}"

rule Star:
    input:
        R1 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R1.fq.gz",
        R2 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R2.fq.gz",
        het_vcf = vcf_path +  "{sample}.snps.het.vcf.gz" 
    output:
        bam = temp(alignment + "{sample}/{sample}_{flowcell}_{lane}.bam")  , 
        bai = temp(alignment + "{sample}/{sample}_{flowcell}_{lane}.bam.bai")
    params:
        rg = " SO=coordinate RGID={flowcell}:{lane} RGLB={sample}.0 RGPL=ILLUMINA RGPU=hiseq RGSM={sample} RGCN=FGCZ ",
        gDir = config['resources']['genomeDir'],
        sjdb = config['resources']['sjdbFile'],
        ref = assembly
    shell:
        """
        
        ## here i kept one combined for STAR not separate rule, just to keep editing it easily
        ## make sure to empty TMPdir usualy snakemake do it after finishing
        
        cd $TMPDIR
        module load jdk 
        
        STAR --runThreadN 15 \
        --twopassMode Basic --genomeDir {params.gDir} --sjdbGTFfile {params.sjdb} --sjdbOverhang 100 \
        --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outSAMmapqUnique 60 --waspOutputMode SAMtag \
        --varVCFfile <(zcat  {input.het_vcf} ) --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random
        
        echo "START adding read groups and sorting"
        #coordinate sorting and add read groups
        
        java -jar picard.jar \
        AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam  O=Aligned_new_sorted.bam {params.rg} 

        ## index 
        samtools index Aligned_new_sorted.bam
        
        ## copy files
        cp $TMPDIR/Aligned_new_sorted.bam {output.bam}
        cp $TMPDIR/Aligned_new_sorted.bam.bai {output.bai}
        
        rm $TMPDIR/Aligned.sortedByCoord.out.bam
        
        """    
        
def list_files(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(sample = wildcards.sample).output.split_fastqs
    all_wildcards = glob_wildcards(os.path.join(checkpoint_output, "{sample}_{flowcell}_{lane}_R1.fq.gz"))
    all_files = []
    for sample, flowcell, lane in zip(all_wildcards.sample, all_wildcards.flowcell, all_wildcards.lane):
        all_files.append(f"{alignment}" + f"{sample}/{sample}_{flowcell}_{lane}.bam")
    return(all_files)


rule sambamba_merge:
    input: list_files
    output: temp(sorted_alignment + "{sample}/{sample}.bam")
    run:
        if len(input) == 1:
            shell("mv {input} {output} \n" + "mv {input}.bai {output}.bai")
        else:
            shell(SAMBAMBA + " merge -t 6 {output} {input}")

rule sambamba_flagstat_1:
    input:
        BAM = sorted_alignment + "{sample}/{sample}.bam"
    output:
        Stats = sorted_alignment + "{sample}/{sample}.stats"
    params:
        threads = "--nthreads 10"
    shell:
        SAMBAMBA + " flagstat -t 3 {params.threads} {input.BAM} > {output.Stats}"

rule mark_duplicates:
    input:
        sorted_alignment + "{sample}/{sample}.bam"
    output:
        bam= dedup_alignment + "{sample}/{sample}.bam",
        metrics= dedup_alignment +"{sample}/{sample}.metrics.txt"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        PICARD_TOOLS + " MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule build_index:
    input:
        dedup_alignment + "{sample}/{sample}.bam"
    output:
        dedup_alignment + "{sample}/{sample}.bam.bai"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        PICARD_TOOLS + " BuildBamIndex I={input} O={output}"

rule picard_metrics:
    input:
        BAM = dedup_alignment + "{sample}/{sample}.bam",
        BAI = dedup_alignment + "{sample}/{sample}.bam.bai",
        ref = assembly
    output:
        metrics = metrics + "{sample}/{sample}.alignment_summary_metrics"
    params:
        prefix = metrics + "{sample}/{sample}"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        LOAD_R +
        PICARD_TOOLS + " CollectMultipleMetrics I={input.BAM} O={params.prefix} R={input.ref}"

rule samtools_idxstats:
    input:
        BAM = dedup_alignment + "{sample}/{sample}.bam",
        BAI = dedup_alignment + "{sample}/{sample}.bam.bai"
    output:
        stats = metrics + "{sample}/{sample}_idxstats"
    shell:
        SAMTOOLS + " idxstats {input.BAM} > {output.stats}"

rule wasp_unique:
    input:
        BAM = dedup_alignment + "{sample}/{sample}.bam",
        BAI = dedup_alignment + "{sample}/{sample}.bam.bai"
    output:
        BAM = wasp_out + "{sample}/{sample}_filter.bam"
    shell:
        """
        
        module load python
        
        python  Filter_wasp_unique.py {input.BAM} {output.BAM}
        

        """
        
rule build_index2:
    input:
        wasp_out + "{sample}/{sample}_filter.bam"
    output:
        wasp_out + "{sample}/{sample}_filter.bam.bai"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        PICARD_TOOLS + " BuildBamIndex I={input} O={output}"
