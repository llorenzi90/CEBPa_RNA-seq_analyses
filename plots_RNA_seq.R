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

require(tidyverse) #esto es una coleccion de packages que ayudan a ordernar la data
require(data.table) #util para leer data rapido y hacer algunas operaciones en tablas 

## ---------------------------
####    Venn diagrams    ####
#existen muchos packages para hacer venn diagrams: https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
#ggvenn es solo uno de tantos...

#si aun no esta instalado correr:
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
#una vez instalado, cargar el paquete:
library(ggvenn)
#Load data
datadir <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/tsv_files/"
datadir <- "/Users/llorenzi/CEBPa/results/DESeq2_results/tsv_files/"
list.files(datadir)#check content of datadir
comparisons <- c("p30vsEV","p42vsEV") 
#read results files into a list
reslist <- lapply(comparisons, function(x)read.table(list.files(datadir,pattern = x,full.names = T),sep = "\t"))
names(reslist) <- comparisons
#take a look at the data:
View(reslist$p30vsEV)
View(reslist$p42vsEV)

#select genes that are deferentially expressed (up and down) for each isoform
cut.padj <- 0.05
#cut.logfc <- #we may also want to filter based on fold change
list.up <- lapply(reslist, function(x)rownames(x)[x$padj<cut.padj&
                                            !is.na(x$padj)& #some genes have padj = NA, remove those (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)  
                                            x$log2FoldChange>0])
list.down <- lapply(reslist, function(x)rownames(x)[x$padj<cut.padj&
                                              !is.na(x$padj)&
                                              x$log2FoldChange<0])
#check how many genes in each isoform
lapply(list.up, length)
lapply(list.down, length)

#define color palette to use in the plots. 
#The one below is a color-blind friendly palette 
cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7") #source http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

#set as working directory the main directory where we want to save the plots
#wdir <- "~/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/"
dir.create("~/tmpdir2")
wdir <- "~/tmpdir2"
setwd(wdir)
#I create a directory for venn diagrams and move to that directory
dir.create("venn_diagrams")
setwd("venn_diagrams/")

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
library(heatmap3) #Sergi suggested to use this package, as it gives nice plots and has a good color palette

#groups of genes
#We want to make a heatmap based on pre-defined groups of genes,
#resulting from the overlap between gene sets defined based on
#differential gene expression analyses 
#between each CEBPA isoform and the empty vector
#we do not want to use the clustering functions of the heatmap package

#So first we have to define our groups of genes.
# I will use the function venn from gplots to calculate the intersection between
#gene sets (this is just a way to do it)
library(gplots)
v.up <- venn(list.up,show.plot = F) #this function also makes a venn diagram (but an ugly one)
                                    #but I indicate not to show the plots and store the results in v.up object
v.down <- venn(list.down,show.plot = F)
#
class(v.up)
"venn"
#this function generates an object of class venn. See more details by running ?venn
#These objects have a list of attributes that we can check like this:
names(attributes(v.up))
#we are intereseted in the "intersections"
names(attributes(v.up)[["intersections"]])

table(attributes(v.up)$intersections$p30vsEV%in%attributes(v.down)$intersections$p42vsEV)
table(attributes(v.down)$intersections$p30vsEV%in%attributes(v.up)$intersections$p42vsEV)
#some genes are upregulated in p30 and downregulated in p42 and vice versa
v.up <- attributes(v.up)$intersections
v.down <- attributes(v.down)$intersections
#define groups:
up_in_both <- v.up$`p30vsEV:p42vsEV`
down_in_both <- v.down$`p30vsEV:p42vsEV`
up_in_p30 <- v.up$p30vsEV
up_in_p30_down_in_p42 <- up_in_p30[up_in_p30%in%v.down$p42vsEV]
up_in_p30 <- up_in_p30[!(up_in_p30%in%up_in_p30_down_in_p42)]
up_in_p42 <- v.up$p42vsEV
up_in_p42_down_in_p30 <- up_in_p42[up_in_p42%in%v.down$p30vsEV]
up_in_p42 <- up_in_p42[!(up_in_p42%in%up_in_p42_down_in_p30)]
down_in_p30 <- v.down$p30vsEV
down_in_p30 <- down_in_p30[!(down_in_p30%in%up_in_p42_down_in_p30)]
down_in_p42 <- v.down$p42vsEV
down_in_p42 <- down_in_p42[!(down_in_p42%in%up_in_p30_down_in_p42)]

#aggregate gene names in a vector 
genes <- c(up_in_both,down_in_both,up_in_p30,down_in_p30,
           up_in_p42,down_in_p42,up_in_p30_down_in_p42,
           up_in_p42_down_in_p30)
length(genes)#5689
length(unique(genes))#5689
length(unique(c(unlist(list.up),unlist(list.down)))) #check that 
#genes in each group are unique and
#that total number matches the original number of unique genes within both lists

#read expression table
#tpm <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.TPM.csv")
tpm <- read.csv("/Users/llorenzi/CEBPa/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.TPM.csv")
#remove LPS samples for now
tpm <- tpm[,grep("LPS",colnames(tpm),invert = T)]
#select genes
table(genes%in%tpm$gene_name)
#dos maneras alternativas de retener los genes de interes:
tpm <- tpm%>%filter(gene_name%in%genes)
tpm <- tpm[tpm$gene_name%in%genes,]

#add group
gl <- list(up_in_both=up_in_both,
           down_in_both=down_in_both,
           down_in_p42=down_in_p42,
           down_in_p30=down_in_p30,
           up_in_p30_down_in_p42=up_in_p30_down_in_p42,
           up_in_p30=up_in_p30,
           up_in_p42_down_in_p30=up_in_p42_down_in_p30,
           up_in_p42=up_in_p42)
tmp <- sapply(gl,function(x)tpm$gene_name%in%x)
group <- apply(tmp,1,function(x)colnames(tmp)[x])
tpm$group <- unlist(group)
tpm$group_index <- match(tpm$group,names(gl))
#tpm <- tpm[order(tpm$group_index),]

tpm$mean_exp <- apply(tpm[,2:10],1,mean)
#sorted_tpm <- tpm %>% group_by(group) %>% arrange( -mean_exp, .by_group = T)
#sorted_tpm_correct_group_order <- sorted_tpm %>%arrange(group_index)
sorted_tpm_correct_group_order <- tpm %>% group_by(group_index) %>% arrange( -mean_exp, .by_group = T)

getwd()
#setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/")
wdir
setwd(wdir)
dir.create("heatmaps")
setwd("heatmaps/")

rowsidecolors <- cbPalette[1:(length(levels(as.factor(sorted_tpm_correct_group_order$group))))][sorted_tpm_correct_group_order$group_index]
pdf("test.pdf")
heatmap3(sorted_tpm_correct_group_order[,2:10],
         Rowv = NA,
         Colv = NA,
         labRow = sorted_tpm_correct_group_order$gene_name,
         RowSideColors = rowsidecolors,
         legendfun=function() showLegend(legend=unique(sorted_tpm_correct_group_order$group),
                                         col=unique(rowsidecolors),cex=1.3),cexCol = 1.5,RowSideLabs = "Group")
dev.off()

