### script to annotate cluster ids

library(argparser, quietly=TRUE)

p <- arg_parser("annotate intron cluster ")
p <- add_argument(p, "--exon", help="exon file with path in proper format chr start end strand geneid genename type")
p <- add_argument(p, "--intron", help="intron cluster id file from top conditional same like exon file format with header")
p <- add_argument(p, "--topfile", help="final file of top conditional output")
p <- add_argument(p, "--outfile", help="output file name ")
p <- add_argument(p, "--output_dir", help="Output directory", default=".")
argv <- parse_args(p)


exon<-read.table(argv$exon) ### 
file<-read.table(argv$intron, header = T) ### other covarites like sex, RIN value etc

colnames(exon)<-c("chr","start","end","strand","gene_id","gene_name","type")
head(exon)

### funtion to annotate

library(dplyr)
library(doParallel)
library(foreach)
library(dplyr)

map_clusters_to_genes <- function(intron_meta, exons_table) {
  gene_df <- foreach (chr=sort(unique(intron_meta$chr)), .combine=rbind) %dopar%  {



   intron_chr <- intron_meta[ intron_meta$chr==chr, ]
    exons_chr <- exons_table[exons_table$chr==chr, ]



   exons_chr$temp <- exons_chr$start
    intron_chr$temp <- intron_chr$end
    three_prime_matches <- inner_join( intron_chr, exons_chr, by=c("temp","strand"))



   exons_chr$temp <- exons_chr$end
    intron_chr$temp <- intron_chr$start
    five_prime_matches <- inner_join( intron_chr, exons_chr, by=c("temp","strand")) 



    all_matches <- rbind(three_prime_matches, five_prime_matches)[ , c("chr.x","start.x", "end.x", "pid", "strand","gene_id","gene_name","type")]
    all_matches <- all_matches[!duplicated(all_matches),]

    all_matches
  }
    
clu_df <- gene_df %>% group_by(pid) %>% summarize(genes=paste(gene_id, collapse = ","), gene_type=paste(type, collapse = ","), gene_id= paste(gene_name, collapse = ","))
class(clu_df) <- "data.frame"
clu_df
    
    }

## annotate intron ids

m <- map_clusters_to_genes(file, exon)

full<-read.table(argv$topfile, header = T) ### 
full$gene<-m$genes[match(full$phe_id , m$pid)]
full$gene_id<-m$gene_id[match(full$phe_id , m$pid)]
full$gene_type<-m$gene_type[match(full$phe_id , m$pid)]

write.table(full,file.path(argv$output_dir,argv$outfile), col.names = T, row.names = F, quote = F, sep = "\t")
