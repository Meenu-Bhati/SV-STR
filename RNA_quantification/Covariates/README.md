## Prepare covariates for eQTL analyses

Run the following command to get peer factors

### 1. Calculate peer factors
```
## script is from Human GTEx consortium https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl

Rscript run_peer.r FINAL_gene_TPM_filter_quan_inv_allgenes.bed.gz prefix 10

```
### 2 Consolidate all covariates and choose the most suitable representative covariates.

```
## script combines all the covariates in one file and further checks and visualize variance explained by each covariates

Rscript Check_and_merge_covariates.r

```

