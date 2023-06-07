#!/bin/bash
module load  gcc/4.8.5 python/3.7.4
snakemake --use-conda --jobs 500 -rp  --latency-wait 120 --rerun-incomplete --snakefile Snakemake_refSTR.py --configfile config.yaml --cluster-config cluster.json  --cluster "bsub -J {cluster.jobname} -n {cluster.ncore} -W {cluster.jobtime} -o {cluster.logi} -R \"rusage[mem={cluster.memo}, scratch={cluster.scratch_mem}]\""
