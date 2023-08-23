## Splicing QTL analyses

### 1. Splicing quantification

To conduct splicing quantification as intron excision ratios, utilize the script "Splicing_quantification.sh." 
Make sure to update the file path accordingly.

### 2. Prepare Covariates

use same [Covariates](https://github.com/Meenu-Bhati/SV-STR/blob/main/RNA_quantification/Covariates/) from eQTL analyses and substitute the expression (TPM) file with the final output from the preceding script (leafcutter output), along with the PCA generated using leafcutter.

### 3. sQTL 

Use [eQTL scripts](https://github.com/Meenu-Bhati/SV-STR/blob/main/eQTL/eQTL_script.sh) from eQTL analyses and substitute the expression (TPM) file with the final output from step 1 (leafcutter output), along with the Covariates generated from previous steps.

Similarily for indepednant SV, STR, SV, Joint SV_STR, and ALL sQTL analyses, Change the VCF file to each respective variant category.
