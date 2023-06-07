## parsing the required files
configfile: "config.yaml"

# Deal with the lsf logfile
import os
os.system("mkdir log_folder")

# Reduced paths
assembly_path =  config['path']['assembly_path']
output_path   =  config['path']['output_path']
assemblies = config['assemblies']
#tools
TRF = config['tools']['trf']
BEDTOOLS = config['tools']['bedtools']

## assemblies 
chrs = list(range(1,30)) + ['X','MT']

rule all:
    input:
        trf_out = expand ( output_path + "/{assembly}/trf_out/{chr}.fa", chr=chrs, assembly=assemblies),
        fix_trf = expand ( output_path + "/{assembly}/fix_trf_out/{chr}.fa" , chr=chrs, assembly=assemblies )
        
rule trf_run:
    input: 
        fasta_file = assembly_path  + "/{assembly}/{chr}.fa"
    output:
        out_fa = output_path + "/{assembly}/trf_out/{chr}.fa"
    params:
        # Default parameters to use with Tandem Repeat Finder
        MATCH_WT = "2",
        MISMATCH_PEN ="7",
        INDEL_PEN ="7",
        P_MATCH ="80",
        P_INDEL ="10",
        MIN_SCORE= "5",
        MAX_PERIOD ="500",
        name = "{chr}"     
    shell:
        TRF + " {input.fasta_file} {params.MATCH_WT} {params.MISMATCH_PEN} {params.INDEL_PEN} {params.P_MATCH} {params.P_INDEL}   {params.MIN_SCORE} {params.MAX_PERIOD} -h -d -l 6 -ngs > {output.out_fa}"

rule fix_trfout:
    input:
        fasta_file = rules.trf_run.output.out_fa
    output:
        out_fix =  output_path + "/{assembly}/fix_trf_out/{chr}.fa"
    conda:
        "py2.yaml"    
    shell:
        """
        python fix_trf_output.py {input.fasta_file} {output.out_fix} 

        """

        
        
        
        
        