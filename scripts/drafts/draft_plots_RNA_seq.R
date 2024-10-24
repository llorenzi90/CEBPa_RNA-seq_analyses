## ---------------------------
##
##
## Purpose of script: Plots for RNA-seq 
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-06-29
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
#Esto es para mantener el formato de los números muy grandes,
#para que, en caso de exportar datos a una tabla R no escriba cosas como '10e9' 
## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ---------------------------
####    Venn diagrams    ####
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
#Load data
datadir <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/"
list.files(datadir)#check content of datadir
comparisons <- c("p30vsEV","p42vsEV") 
#read results files into a list
reslist <- lapply(comparisons, function(x)read.csv(list.files(datadir,pattern = x,full.names = T),sep = ";"))
names(reslist) <- comparisons
View(reslist$p30vsEV)
View(reslist$p42vsEV)
#select genes that are differentially expressed (up and down) for each isoform
cut.padj <- 0.05
#cut.logfc <- 
list.up <- lapply(reslist, function(x)x$X[x$padj<cut.padj&
                                            !is.na(x$padj)&
                                            x$log2FoldChange>0])
list.down <- lapply(reslist, function(x)x$X[x$padj<cut.padj&
                                              !is.na(x$padj)&
                                              x$log2FoldChange<0])
lapply(list.up, length)
lapply(list.down, length)

#colors
cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7") #source http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
pdf("upgenes.pdf")
ggvenn(
  list.up, 
  fill_color = cbPalette[1:2], #change this to set other colors 
  #fill_color = c( "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4,
) + ggtitle("Upregulated genes")
dev.off()

pdf("downgenes.pdf")
ggvenn(
  list.down, 
  fill_color = cbPalette[1:2], #change this to set other colors 
  #fill_color = c( "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4,
) + ggtitle("Downregulated genes")
dev.off()


####      Heat maps     ####
install.packages("heatmap3")
library(heatmap3)

#groups of genes
library(gplots)
v.up <- venn(list.up)
v.down <- venn(list.down)
table(attributes(v.up)$intersections$p30vsEV%in%attributes(v.down)$intersections$p42vsEV)
#some genes are upregulated in p30 and downregulated in p42 and vice versa
table(attributes(v.down)$intersections$p30vsEV%in%attributes(v.up)$intersections$p42vsEV)

up_in_both <- attributes(v.up)$intersections$`p30vsEV:p42vsEV`
down_in_both <- attributes(v.down)$intersections$`p30vsEV:p42vsEV`
up_in_p30 <- attributes(v.up)$intersections$p30vsEV
up_in_p30_down_in_p42 <- up_in_p30[up_in_p30%in%(attributes(v.down)$intersections$p42vsEV)]
up_in_p30 <- up_in_p30[!(up_in_p30%in%up_in_p30_down_in_p42)]
up_in_p42 <- attributes(v.up)$intersections$p42vsEV
up_in_p42_down_in_p30 <- up_in_p42[up_in_p42%in%attributes(v.down)$intersections$p30vsEV]
up_in_p42 <- up_in_p42[!(up_in_p42%in%up_in_p42_down_in_p30)]
down_in_p30 <- attributes(v.down)$intersections$p30vsEV
down_in_p30 <- down_in_p30[!(down_in_p30%in%up_in_p42_down_in_p30)]
down_in_p42 <- attributes(v.down)$intersections$p42vsEV
down_in_p42 <- down_in_p42[!(down_in_p42%in%up_in_p30_down_in_p42)]
length(unique(c(up_in_both,down_in_both,up_in_p30,down_in_p30,
                up_in_p42,down_in_p42,up_in_p30_down_in_p42,
                up_in_p42_down_in_p30)))
length(unique(c(unlist(list.up),unlist(list.down))))
genes <- unique(c(unlist(list.up),unlist(list.down)))
#read expression table
library(data.table)
tpm <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.TPM.csv")
#remove LPS samples for now
tpm <- tpm[,grep("LPS",colnames(tpm),invert = T)]
#select genes
tpm <- tpm%>%filter(gene_name%in%genes)
table(genes%in%tpm$gene_name)
#FALSE  TRUE 
#2  5687
#two genes are not in the count table,
#let's check what these genes are
genes[!genes%in%tpm$gene_name]
#[1] "ago-02" "ago-01"
#It's Excel's fault! It converts some gene names to dates...
#Note: never save files in Excel... preferably never use Excel...
#I will run again DESeq2 and save as tsv