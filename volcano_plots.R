## ---------------------------
##
##
## Purpose of script:volcano plots CEBPa RNA-seq
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-04-04
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
library(ggrepel)
## ---------------------------
padj.cutoff <- 0.05
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/tsv_files/")
ress=list.files(pattern = "DESeq_results")
#palette of colors to use
cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7") #source http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/volcano_plots/"
padjco=0.05
Ntop=5

for (resf in ress) {
  da=read.table(resf)
  nam=gsub("_DESeq_results.txt","",resf)
  da <- da[!is.na(da$padj),]
  
  da$diffexpressed <- ifelse(da$padj<=0.05,ifelse(da$log2FoldChange>0,"UP","DOWN"),"NO")
  #da$diffexpressed <- ifelse(da$padj<=0.05&abs(da$log2FoldChange)>=0.58,
  #ifelse(da$log2FoldChange>0,"UP","DOWN"),"NO") #this were the paraneters for old_plots
  print(table(da$diffexpressed))
  
  #Select genes to show and convert padj=0 to next min value
  da$genes_names_to_show <- rownames(da)
  
  da$genes_names_to_show[order(da$padj)][6:nrow(da)] <- NA
  table(da$padj==0)
  da$padj[da$padj==0] <- min(da$padj[da$padj!=0])
  
  
  
  #base ggplot
  p <- ggplot(data=da, 
              aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed)) + 
    geom_point() +
    theme_classic()
  #use cb palette:
  pf=p + geom_text_repel(data = da ,
                      aes(label=genes_names_to_show),
                      max.overlaps = 25) +   
    scale_color_manual(values = c(cbPalette[7],cbPalette[1],cbPalette[6]))+
    ggtitle(nam)+
    theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))
  pdf(paste0(outdir,nam,"_volcano.pdf"))
  print(pf)
  dev.off()  
  tiff(paste0(outdir,nam,"_volcano.tiff"),units="in", width=8, height=6, res=300, compression = 'lzw')
  print(pf)
  dev.off()  
  }

