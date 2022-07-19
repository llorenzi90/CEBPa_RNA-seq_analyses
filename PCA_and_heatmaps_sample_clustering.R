options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
library(DESeq2)
library(ggplot2)
countdata=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.txt",header = T)
rownames(countdata) <- countdata$gene_name
countdata <- countdata[,-1]
#quick PCA
libsizes <- colSums(countdata)
scaled_counts <- apply(countdata,1, function(x)x/libsizes*1000000)
scaled_counts <- scaled_counts[,colVars(scaled_counts)!=0]
pca <- prcomp(scaled_counts,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
spt=strsplit(dtp$samples,split = "\\.")
dtp$vector=sapply(spt,function(x)x[1])
#dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=sapply(spt,function(x)any(grepl("LPS",x)))
dtp$rep=sapply(spt,function(x)ifelse(x[2]=="LPS",return(x[3]),return(x[2])))
lps=dtp$LPS
dtp$LPS[!lps] <- "UT"
dtp$LPS[lps] <- "LPS"

summary_alignments=read.csv('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/QC/summary_alignments.csv')
summary_alignments$name <- gsub("-",".",summary_alignments$name)
colnames(summary_alignments)[1:4] <- c("SampleID","Condition","Treatment","Replicate")
table(summary_alignments$SampleID%in%dtp$samples)
summary_alignments$Treatment <- dtp$LPS[match(summary_alignments$SampleID,
                                              dtp$samples)]
summary_alignments$pairs_mapped_to_genes <- libsizes[match(summary_alignments$SampleID,
                                                           dtp$samples)]
summary_alignments$SampleID=paste(summary_alignments$Condition,summary_alignments$Treatment,
                                  summary_alignments$Replicate,sep = "_")
write.table(summary_alignments,'/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/QC/summary_alignments.reformatted.csv',row.names = F)


#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/")
dir.create("PCA_plots")

pdf("PCA_plots/PCA_CPM.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()
pdf("PCA_plots/PCA_CPM.colLPS.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector),size=4) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))+theme(text = element_text(size = 20))
dev.off()

##Add pca using 500 most variable genes (this is what DESeq2 package does by default and
#what I have been doing for other analyses) as a comparison
rv <- rowVars(t(scaled_counts))
ntop=500 #select only the top 500 genes with highest variance
#This is optional, is the default of the DESeq2 package
select <- order(rv, decreasing = TRUE)[1:ntop]
top500variable <- scaled_counts[,select]
pca <- prcomp(top500variable,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
spt=strsplit(dtp$samples,split = "\\.")
dtp$vector=sapply(spt,function(x)x[1])
#dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=sapply(spt,function(x)any(grepl("LPS",x)))
dtp$rep=sapply(spt,function(x)ifelse(x[2]=="LPS",return(x[3]),return(x[2])))
lps=dtp$LPS
dtp$LPS[!lps] <- "UT"
dtp$LPS[lps] <- "LPS"
#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

pdf("PCA_plots/PCA_CPM.top500moreVarGenes.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()
pdf("PCA_plots/PCA_CPM.top500moreVarGenes.colLPS.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector),size=4) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))+theme(text = element_text(size = 20))
dev.off()

#####log2

pca <- prcomp(log2(top500variable+0.01),scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
spt=strsplit(dtp$samples,split = "\\.")
dtp$vector=sapply(spt,function(x)x[1])
#dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=sapply(spt,function(x)any(grepl("LPS",x)))
dtp$rep=sapply(spt,function(x)ifelse(x[2]=="LPS",return(x[3]),return(x[2])))
lps=dtp$LPS
dtp$LPS[!lps] <- "UT"
dtp$LPS[lps] <- "LPS"
#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

pdf("PCA_plots/PCA_CPM.top500moreVarGenes.log2+0.01.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()
pdf("PCA_plots/PCA_CPM.top500moreVarGenes.colLPS.log2+0.01.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector),size=4) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))+theme(text = element_text(size = 20))
dev.off()

######heatmap
sampleDists <- dist(scaled_counts)
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("cluster_heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf("cluster_heatmaps/pheatmap_sample_clustering_allgenes.pdf")

print(p)
dev.off()



#log2
sampleDists <- dist(log2(scaled_counts+0.01))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 12)
pdf("cluster_heatmaps/pheatmap_sample_clustering_allgenes.log2+0.01.pdf")

print(p)
dev.off()

#top500variable

sampleDists <- dist(top500variable)
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("cluster_heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf("cluster_heatmaps/pheatmap_sample_clustering_top500variable.pdf")

print(p)
dev.off()



#log2
sampleDists <- dist(log2(top500variable+0.01))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 12)
pdf("cluster_heatmaps/pheatmap_sample_clustering_top500variable.log2+0.01.pdf")

print(p)
dev.off()
