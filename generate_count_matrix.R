setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/htseq_count/data")
library(data.table)
countfiles=list.files(pattern = ".htseq_count" )
htseq_counts_all=read.table(countfiles[1])
for (fi in countfiles[-1]) {
  htseq_counts_all <- cbind(htseq_counts_all,read.table(fi))
}

non_mapped_stats <- htseq_counts_all[grep("ENSMUSG",rownames(htseq_counts_all),invert = T),]
htseq_counts <-  htseq_counts_all[grep("ENSMUSG",rownames(htseq_counts_all)),]
total_sequenced_reads <- apply(htseq_counts_all,2,sum)
total_gene_counts <- apply(htseq_counts,2,sum)
total_gene_counts/total_sequenced_reads*100
total_non_mapped_counts=apply(non_mapped_stats,2, sum)
all_stats=rbind("__total_non_mapped_counts"=total_non_mapped_counts,non_mapped_stats)
all_stats=rbind("__total_gene_counts"=total_gene_counts,all_stats)
all_stats=rbind("__total_sequenced_counts"=total_sequenced_reads,all_stats)
write.csv(all_stats,"../htseq_count_stats.csv")
#write.csv(cbind(gene_id=rownames(htseq_counts),htseq_counts),"../../expression_files/all_samples_gencodevM27.counts.csv",quote = F,row.names = F)
