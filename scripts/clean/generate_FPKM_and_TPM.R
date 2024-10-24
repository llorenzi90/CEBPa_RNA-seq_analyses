## ---------------------------
##
##
## Purpose of script: generate FPKM/TPM expression files 
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

#require(tidyverse)
#require(data.table)

## ---------------------------

countdata <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv")

### Normalize data ###
#remove rownames to work for now
gene_names <- countdata$gene_name
countdata <- countdata[,-1]
 
#read gene lengths and gene ids/gene names table
gene_lengths <- read.table("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/references/mouse/gencode/gencode.vM27.gene_lengths.tsv")
geneid_gene_names <- read.table("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/references/mouse/gencodeVM27_annot_ENSEMBLgeneID_to_unduplicated_gene_name.tsv",header = T)
table(gene_lengths$V1==geneid_gene_names$gene_id)
table(geneid_gene_names$gene_name==gene_names)
#check if the file is already in the right order 
# 
##TPM
# 
##FPKM
#These three metrics attempt to normalize for sequencing depth and gene length. Here’s how you do it for RPKM (or FPKM):
#   
#1)Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
pmsf <- apply(countdata, 2,sum)/1000000
#2)Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
rpm <- t(apply(countdata,1,function(x)x/pmsf))
head(rpm)
#3)Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
kilobases <- gene_lengths$V2/1000
FPKM <- apply(rpm,2,function(x)x/kilobases )
# 
# ########compare to result using fpkm function from DESeq2#######
# library(DESeq2)
# #extrafont::loadfonts()
# countdata=read.csv("data/all_samples_gencodevM27.counts.csv")
# rownames(countdata) <- countdata$gene_id
# countdata <- countdata[,-1]
# name=colnames(countdata)
# tr=strsplit(name,split = "\\.")
# vector=sapply(tr, function(x)return(x[1]))
# rep=sapply(tr, function(x)return(x[length(x)]))
# LPS=sapply(tr, function(x)return(any(grepl("LPS",x))))
# coldata=data.frame(name,vector,LPS,rep)
# rownames(coldata) <- coldata$name
# table(rownames(coldata)==colnames(countdata))
# coldata$vector <- as.factor(coldata$vector)
# coldata$vector
# class(coldata$LPS)
# coldata$LPS <- as.factor(coldata$LPS)
# coldata$LPS
# 
# ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
#                                  colData = coldata,
#                                  design = ~0+ vector + LPS + vector:LPS)
# 
# mcols(ddsMat)$basepairs <- gene_lengths$V2
# test <- fpkm(object=ddsMat,robust = F)
# 
FPKM <- cbind(gene_name=rownames(FPKM),FPKM)
write.csv(FPKM,"/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.FPKM.csv",row.names = F,quote = F)
# 
# 
#####TPM####
countdata <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv")
rownames(countdata) <- countdata$gene_name
countdata <- countdata[,-1]

#1)Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
rpk <- apply(countdata, 2, function(x)x/kilobases)
#2)Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf.1 <- apply(rpk, 2, sum)/1000000
#3)Divide the RPK values by the “per million” scaling factor.
TPM <- t(apply(rpk, 1, function(x)x/pmsf.1))
# 
TPM <- cbind(gene_name=rownames(TPM), TPM)
write.csv(TPM,"/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.TPM.csv",row.names = F,quote = F)
