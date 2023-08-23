## Prepare covariates for eQTL analyses

Run follwing to commond to get peer 

### 1. calculate peer factors
```
## script is from Human GTEx consortium https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl

Rscript run_peer.r FINAL_gene_TPM_filter_quan_inv_allgenes.bed.gz prefix 10

```
### 2 Consolidate all covariates and choose the most suitable representative covariates.

```
## script combine all the covriates in one file and further check and visualize variance explained by each covarites

Rscript Check_and_merge_covriates.r

```

