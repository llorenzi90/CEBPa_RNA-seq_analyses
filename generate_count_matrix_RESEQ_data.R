#setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/htseq_count/data")
setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/RESEQ_data/data/htseq_count/")
library(data.table)
countfiles=list.files(pattern = "out" )
htseq_counts_all=read.table(countfiles[1],header = T,row.names = 1)
for (fi in countfiles[-1]) {
  htseq_counts_all <- cbind(htseq_counts_all,read.table(fi,header = T,row.names = 1))
}
library(tidyverse)
non_mapped_stats <- htseq_counts_all[grep("ENSMUSG",rownames(htseq_counts_all),invert = T),]
htseq_counts <-  htseq_counts_all[grep("ENSMUSG",rownames(htseq_counts_all)),]
total_sequenced_reads <- apply(htseq_counts_all,2,sum)
total_gene_counts <- apply(htseq_counts,2,sum)
total_gene_counts/total_sequenced_reads*100
total_non_mapped_counts=apply(non_mapped_stats,2, sum)
all_stats=rbind("__total_non_mapped_counts"=total_non_mapped_counts,non_mapped_stats)
all_stats=rbind("__total_gene_counts"=total_gene_counts,all_stats)
all_stats=rbind("__total_sequenced_counts"=total_sequenced_reads,all_stats)
#write.csv(all_stats,"htseq_count_stats.csv")
#write.csv(cbind(gene_id=rownames(htseq_counts),htseq_counts),"../expression_files/all_samples_gencodevM27.counts.csv",quote = F,row.names = F)

all_stats_previous_seq <- read.csv("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/htseq_count/htseq_count_stats.csv",row.names = 1)
all_stats/all_stats_previous_seq
View(all_stats/all_stats_previous_seq)


#A lot of aligned reads are not assigned to any gene!! 
#check correlation between reseq data and previous data
old_expression <- read.csv("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/htseq_count/data/all_samples_gencodevM27.counts.csv",row.names = 1)
colnames(old_expression) <- paste0(colnames(old_expression),"_prev")

all_data <- cbind(htseq_counts,old_expression)
library(DESeq2)
vsd <- vst(as.matrix(all_data),blind = T)

sampleDists <- dist(t(vsd))
sampleDists

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#an alternative way to calculate distance:
library("PoiClaClu")
poisd <- PoissonDistance(t(vsd))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste(colnames(all_data))
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plot(vsd[,"p30.LPS.R1"],vsd[,"p30.LPS.R1_prev"])
plot(vsd[,"p30.R1"],vsd[,"p30.R1_prev"])
plot(vsd[,"p30.R1"],vsd[,"p30.LPS.R1"])
plot(vsd[,"p30.R1_prev"],vsd[,"p30.LPS.R1_prev"])

#there are some genes that are more expressed in the new sequencing
#what are those genes?
rat <- vsd[,"p30.LPS.R1"]/vsd[,"p42.LPS.R1_prev"]
fil <- rat>=1.2
plot(vsd[,"p30.LPS.R1"],vsd[,"p30.LPS.R1_prev"])
points(vsd[fil,"p30.LPS.R1"],vsd[fil,"p30.LPS.R1_prev"],col="red")

#are the same genes up also in other samples?
plot(vsd[,"p42.R1"],vsd[,"p42.R1_prev"])
points(vsd[fil,"p42.R1"],vsd[fil,"p42.R1_prev"],col="red")

#what are those genes?
up_genes_newseq <- rownames(vsd)[fil]

nofeature_reads <- read.table("~/nofeature.sam",nrows = 500000,sep = "\t",quote = "",fill = T)
sort(table(nofeature_reads$V2),decreasing = T)
#select: proper pairs (147 ,    99,    163  , 83), mates in fwd strand (99 and 163)
nofeature_reads_f <- nofeature_reads %>%filter(V2%in%c(99,163))
#strand 99 = "-", 163 = "+"
#generate bed file to display in IGV:
nofeature_reads_f$strand="+"
nofeature_reads_f$strand[nofeature_reads_f$V2=="99"] <- "-"

head(nofeature_reads_f)
nofeature_reads_f$end <- nofeature_reads_f$V4 + abs(nofeature_reads_f$V9)
nofeature_reads_bed <- nofeature_reads_f[,c(3,4,26,1,5,25)]
write.table(nofeature_reads_bed,"~/test.nofeature_reads.bed",quote = F,row.names = F,col.names = F,sep = "\t")

#some regions to check
#chr12:100,168,171-100,178,201

#
ENSMUSG00000062604_sam=fread("~/ENSMUSG00000062604.sam",sep = "\t",fill = T)
table(ENSMUSG00000062604_sam$V2)
ENSMUSG00000062604_sam_f <- ENSMUSG00000062604_sam %>%filter(V2%in%c(99,163))
ENSMUSG00000062604_sam_f$strand="+"
ENSMUSG00000062604_sam_f$strand[ENSMUSG00000062604_sam_f$V2=="99"] <- "-"

head(ENSMUSG00000062604_sam_f)
ENSMUSG00000062604_sam_f$end <- ENSMUSG00000062604_sam_f$V4 + abs(ENSMUSG00000062604_sam_f$V9)
ENSMUSG00000062604_bed <- ENSMUSG00000062604_sam_f[,c(3,4,26,1,5,25)]
write.table(ENSMUSG00000062604_bed,"~/ENSMUSG00000062604.bed",quote = F,row.names = F,col.names = F,sep = "\t")
