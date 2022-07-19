summary_alignments.reformatted <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/QC/summary_alignments.reformatted.csv", sep="")
View(summary_alignments.reformatted)

library(ggplot2)

library(tidyr)
gd=summary_alignments.reformatted %>% tidyr::gather(key = "step",value = "Nreads",c(5,7,12))
gdlevels(as.factor(gd$step))
g=ggplot(gd,aes(x=factor(step,
                          levels=c("total_sequenced_read_pairs",
                                   "properly_aligned_pairs",
                                   "pairs_mapped_to_genes")),y=Nreads/1000000))+
 geom_boxplot()+theme_classic()+
  ylab("M paired reads")+ 
  xlab("")+
  scale_x_discrete(labels=c("sequenced","properly aligned","mapped to genes"))+
  theme(text = element_text(size = 20))

pdf("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/QC/sequencing_numbers_boxplot.pdf")  
print(g)
dev.off()

tiff("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/QC/sequencing_numbers_boxplot.tiff",res = 300,width = 2000,height = 2000)  
print(g)
dev.off()
