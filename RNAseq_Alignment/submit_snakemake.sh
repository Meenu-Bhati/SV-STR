#!/bin/bash
module load gcc/6.3.0 ### we are using this python and gcc for  gdc-fastq-splitter/myenv/bin/gdc-fastq-splitter"
module load python/3.8.5

snakemake --jobs 5100 -rp --latency-wait 40 --keep-going --rerun-incomplete --snakefile Snakemake.py --cluster-config cluster.json --cluster "bsub -J {cluster.jobname} -n {cluster.ncore} -W {cluster.jobtime} -oo {cluster.logi} -R \"rusage[mem={cluster.memo},scratch={cluster.scratch_mem}]\""
