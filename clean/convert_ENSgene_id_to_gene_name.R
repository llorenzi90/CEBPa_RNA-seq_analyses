## ---------------------------
##
##
## Purpose of script: convert ENSEMBL gene ids to gene names 
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-05-28
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ---------------------------

countdata <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/htseq_count/data/all_samples_gencodevM27.counts.csv")
#1st translate gene ids to gene names
library(rtracklayer)
geneannot <- readGFF("~/share/Cuartero Group/CUARTERO GROUP/references/mouse/gencode/gencode.vM27.annotation.gtf")

table(countdata$gene_id%in%geneannot$gene_id)

geneannot_genes <- geneannot%>%filter(type=="gene")

length(unique(geneannot_genes$gene_id))
#gene ids are not duplicated, this is as expected
length(unique(geneannot_genes$gene_name))
#some gene names are duplicated, so I can't assign them as rownames
#I will add numbering (i.e GENEA_1, GENEA_2, GENEA_3...)
#to duplicated gene names 
#This is just to keep track of possible weird things, but 
#most likely these genes with duplicated gene names are not
# relevant for our analyses, and they likely have 
#aprox 0 counts as they are genes duplicated in different regions

####deduplication of gene names:####
#select a matrix with gene id and corresponding gene name
geneid_gene_names <-  geneannot_genes[match(countdata$gene_id,
                                            geneannot_genes$gene_id),c("gene_id","gene_name")]
#subselect those genes with repeated gene names
geneid_gene_names_dup <- geneid_gene_names%>%filter(duplicated(gene_name))

#keep track of the order of original index when 
#gene names are alphabeticaly sorted (this is because next 
#function will sort names by alphabetic order)
ord <- order(geneid_gene_names_dup$gene_name)
head(geneid_gene_names_dup[ord,])

#count how many times each duplicated gene is duplicated 
ta <- table(geneid_gene_names_dup$gene_name)
table(ta)
# ta
# 1   2   3   4   5   6   7   8   9  11  12  58  67 130 329 337 485 
# 83   6   2   1   1   2   1   1   1   3   1   1   1   1   1   1   1 
#most genes are duplicated only one time,
#but some genes are duplicated hundreds of times

#check that the gene names in original matrix when
#sorted with the ord index are the same than the names in ta
table(unique(geneid_gene_names_dup$gene_name[ord])==names(ta))

#generate gene suffixes for each gene and the corresponding 
#modified gene names
ta_vect <- sapply(ta, function(x)return(paste0("_",seq(1:x))))
modified_names <- paste0(geneid_gene_names_dup$gene_name[ord],unlist(ta_vect))

#assign the modified gene names to the original matrix in the 
#correct order
geneid_gene_names_dup[ord,"modified_gene_names"] <- modified_names

#now add those modified names to the whole matrix in places where duplicated
#genes are found:
geneid_gene_names$gene_name[duplicated(geneid_gene_names$gene_name)] <- geneid_gene_names_dup$modified_gene_names 
#####

# assign the corrected gene names as rownames for the countdata
rownames(countdata) <- geneid_gene_names$gene_name[match(countdata$gene_id,
                                                         geneid_gene_names$gene_id)]


countdata$gene_id <- rownames(countdata)
colnames(countdata)[1] <- "gene_name"
write.csv(countdata,"/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv",row.names = F,quote = F)
write.table(countdata,"/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.txt",row.names = F,quote = F,sep = "\t")
write.table(geneid_gene_names,"/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/references/mouse/gencodeVM27_annot_ENSEMBLgeneID_to_unduplicated_gene_name.tsv",quote = F,row.names = F)
