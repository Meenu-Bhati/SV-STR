file<-read.csv("FINAL_gene_TPM.tsv", header=T, sep="\t")


##
### create PCA
library(stats)
library(ggplot2)
install.packages("ggfortify")
library(ggfortify)

#transpose the matrix 
options(repr.plot.width=8, repr.plot.height=5)
M <- t(file[,c(7:82)]) ## select sample coloums with TPM values 
# transform the counts to log2 scale 
T <- log2(M + 1)
#compute PCA 
pcaResults_log <- prcomp(T)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
autoplot(pcaResults_log,label = TRUE, label.size = 3,  shape = FALSE) + theme_bw(base_size = 12)
autoplot(pcaResults_log) + theme_bw(base_size = 12)

### subset the data based on expression outlier if any
remove<-c("sample50")
file_sampleout = file[,!(names(file) %in% remove)]


### remove gene with expression lowers than 0.2 TPM in 20% samples

file_updated <- file_sampleout [rowSums(file_sampleout[,c(7:81)] <= 0.2) <= 16, ]
nrow(file_updated)


## quantile standerization per animal
file_filter <- as.matrix(file_updated[,c(7:81)])
library(RColorBrewer)
n <- 10000
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
plot(density(file_filter[,1]),xlim = c(-1, 1000), col = "red")
#for(i in 2:40) lines(density(file_filter[,i ],xlim = c(-1, 1000)), col = col_vector[i])  #

### per gene normal distribution
plot(density(file_filter[1,]), col = "red" )
for(i in 2:10000) lines(density(file_filter[i, ]), col = col_vector[i])

## this shows that data is not normalized


## quantile normalise data
library(preprocessCore)

file_filter <- as.matrix(file_updated[,c(7:81)])  ## number of coloums/samples with expression 
file_filter_norm <- normalize.quantiles(file_filter, copy = TRUE)
head(file_filter )
nrow(file_filter_norm)

rownames(file_filter_norm)<-rownames(file_filter)
colnames(file_filter_norm)<-colnames(file_filter)
head(as.data.frame(file_filter_norm))

plot(density(file_filter_norm[,1]), col = "red")


## Best post https://bsunfullsite.github.io/post/quantile-nomrlization-and-inverse-normal-transform/quantile-normalization-and-inverse-normal-transform/

library(RNOmni)

library(RColorBrewer)
n <- 10000

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
plot(density(file_filter_norm[,1]),xlim = c(-1, 1000), col = "red")
for(i in 2:40) lines(density(file_filter_norm[,i ]), col = col_vector[i])  #


## Best post https://bsunfullsite.github.io/post/quantile-nomrlization-and-inverse-normal-transform/quantile-normalization-and-inverse-normal-transform/
expr.int = t(apply(file_filter_norm, 1, RankNorm ))

head(expr.int)
library(RColorBrewer)
n <- 60

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

### per gene normal distribution
plot(density(expr.int[1,]), col = "red" )
for(i in 2:40) lines(density(expr.int[i, ]), col = col_vector[i])


## add info file
file_quan_int <- merge(file_updated[,c(1:6)], expr.int , by="row.names", all=TRUE) 
file_quan_int<- file_quan_int[,-c(1)]

file_quan_int_order<- file_quan_int[ order( file_quan_int[,1], file_quan_int[,2]),]

head(file_quan_int_order)

write.table(file_quan_int_order, "FINAL_gene_TPM_filter_quan_inv_allgenes.bed", col.names = TRUE, quote=F, sep="\t", row.names=F) 

