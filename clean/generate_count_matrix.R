setwd("htseq_count") #replace with dir where htseq_count output files are located

countfiles=list.files(pattern = ".htseq_count" )
samples=gsub(".htseq_count","",countfiles)                                                                                          

htseq_counts_all=read.table(countfiles[1],row.names = 1)

colnames(htseq_counts_all) <- samples[1]

for (i in 2:length(countfiles)) {
  tmp=read.table(countfiles[i],row.names = 1)
  
  colnames(tmp) <- samples[i]
  htseq_counts_all <- cbind(htseq_counts_all,tmp)
}
non_mapped_stats <- htseq_counts_all[grep("ENSM",rownames(htseq_counts_all),invert = T),]
htseq_counts <-  htseq_counts_all[grep("ENSM",rownames(htseq_counts_all)),]
total_sequenced_reads <- apply(htseq_counts_all,2,sum)
total_gene_counts <- apply(htseq_counts,2,sum)
total_gene_counts/total_sequenced_reads*100
total_non_mapped_counts=apply(non_mapped_stats,2, sum)
all_stats=rbind("__total_non_mapped_counts"=total_non_mapped_counts,non_mapped_stats)
all_stats=rbind("__total_gene_counts"=total_gene_counts,all_stats)
all_stats=rbind("__total_sequenced_counts"=total_sequenced_reads,all_stats)
write.csv(cbind(gene_id=rownames(htseq_counts),htseq_counts),"all_samples_gencodevM27.counts.csv",quote = F,row.names = F)
write.csv(t(all_stats),"htseq_stats.csv")
