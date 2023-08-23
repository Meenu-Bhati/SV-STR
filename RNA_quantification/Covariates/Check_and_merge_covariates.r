library(variancePartition)

expr<-read.table(gzfile("FINAL_gene_TPM_filter_quan_inv_allgenes.bed.gz"), header = T, comment.char = "")


rownames(expr)<-expr$gene

expr<-expr[,-c(1:6)] ## just keep sample coloums

cov<-read.csv("testis_cov.csv") ## such as RIN values etc
PCA<-read.table("All_chr_75_LDpruned_pca.eigenvec") ## generate via plink
colnames(PCA)<-c("fam","ID",paste0(rep("PC", 3), rep(1:20, each = 1)))


### merge other cofactors and PCA
all_data<-merge(PCA[,-c(1)], cov, by.x = "ID", by.y = "ID") ## First coloum is family id which is not used


library(tidyverse)

peer<-read.table("FINAL_gene_TPM_filter_quan_inv_allgenes.PEER_covariates.txt", header=T) ## generated with run peer script
peer$ID<-as.factor(peer$ID)
trans_peer<-t(peer)

names(trans_peer)<-trans_peer[1,]

peer_correct<-as.data.frame(trans_peer)

names(peer_correct)<- peer_correct[1,]
colnames(peer_correct)<-paste0(rep("peer", 3), rep(1:10, each = 1))
peer_correct <- peer_correct[-1,]
peer_correct$ID <- rownames(peer_correct)

peer_correct[,c(1:10)]<-lapply(peer_correct[c(1:9)], function(x) if(is.character(x)) as.numeric(x) else x)    

## finallly combine all covariates

all_data_peer<-merge(all_data, peer_correct, by.x = "ID", by.y = "ID")

rownames(all_data_peer)<- all_data_peer$ID

                               
### check varaince explained by each covariates per gene
## you can add mutiple covariates in the form based on your data and chcek each step below to choose best fit covariates
                               
form <- ~ peer1 + peer2 + peer3 + peer4 + peer5 + peer6  + PC1 + PC2 + PC3 +  RIN + age 

varPart <- fitExtractVarPartModel( as.matrix(expr), form, all_data_peer )
vp <- sortCols( varPart )
plotVarPart( vp )
fig <- plotVarPart( vp )
ggsave(file, fig)

### chcek Collinearity
                               
res <-fitVarPartModel(as.matrix(expr), form, all_data_peer )
colinearityScore( res[[1]] )

# evaluate the collinearity score on the first model fit
# this reports the correlation matrix between coefficient estimates
# for fixed effects
# the collinearity score is the maximum absolute correlation value
# If the collinearity score > .99 then the variance partition
# estimates may be problematic
# In that case, a least one variable should be omitted
colinearityScore( res[[1]] )
                               

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C = canCorPairs( form, all_data_peer)
# Plot correlation matrix
plotCorrMatrix( C )

## extract residuals 

form <- ~ peer1 + peer2 + peer3 + peer4 + peer5 + PC1 + PC2 + PC3 + RIN + age 
varPart <- fitExtractVarPartModel( as.matrix(expr), form, all_data_peer )
vp <- sortCols( varPart )
plotVarPart( vp )
fig <- plotVarPart( vp )
file="variance_partition.pdf"
ggsave(file, fig)
                           
