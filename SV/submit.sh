#!/bin/bash
## gcc version gcc/4.8.5 snakemake version 6.10.0
#module load python/3.6.5
snakemake --use-conda --jobs 500 -rp --latency-wait 120 --rerun-incomplete --snakefile Snakemake_new.py --configfile config.yaml --cluster-config cluster.json  --cluster "bsub -J {cluster.jobname} -n {cluster.ncore} -W {cluster.jobtime} -o {cluster.logi} -R \"rusage[mem={cluster.memo}, scratch={cluster.scratch_mem}]\""
