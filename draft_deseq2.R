## ---------------------------
##
##
## Purpose of script: DESeq2 DGEA of CBPa isoforms upon induction with LPS  
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-05-19
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
setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/htseq_count")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

library(DESeq2)
#extrafont::loadfonts()
countdata=read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv")
rownames(countdata) <- countdata$gene_id
countdata <- countdata[,-1]
name=colnames(countdata)
tr=strsplit(name,split = "\\.")
vector=sapply(tr, function(x)return(x[1]))
rep=sapply(tr, function(x)return(x[length(x)]))
LPS=sapply(tr, function(x)return(any(grepl("LPS",x))))
coldata=data.frame(name,vector,LPS,rep)
rownames(coldata) <- coldata$name
table(rownames(coldata)==colnames(countdata))
coldata$vector <- as.factor(coldata$vector)
coldata$vector
class(coldata$LPS)
coldata$LPS <- as.factor(coldata$LPS)
coldata$LPS

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~0+ vector + LPS + vector:LPS)

#filter non-informative rows
#removing rows of the DESeqDataSet that have no counts, or only a single count across all samples. Additional weighting/filtering to improve power is applied at a later step in the workflow.
nrow(ddsMat)
## [1] 58294
keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]
nrow(ddsMat)
## [1] 31604
ddsMat <- estimateSizeFactors(ddsMat)

#Data visualization
#transform with vst
vsd<- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)

#blind = FALSE, means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).

#compare log2 transformation with vst
library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = EV.LPS.R1, y = EV.LPS.R2)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#calculate sample-sample distances
sampleDists <- dist(t(assay(vsd)))
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
poisd <- PoissonDistance(t(counts(ddsMat)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( ddsMat$name)
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#PCA
plotPCA(vsd,intgroup=c("vector","LPS"))

###########################differential expression########################
#the design was already indicated in 
#ddsMat <- DESeqDataSetFromMatrix(countData = countdata,colData = coldata,
#design = ~0+ vector + LPS + vector:LPS)

#ExploreModelMatrix is a cool tool to help to generate the contrasts of interest:

#BiocManager::install("ExploreModelMatrix")
library("ExploreModelMatrix")

vd <- VisualizeDesign(sampleData = coldata, 
                      designFormula = ~0+ vector + LPS + vector:LPS, 
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
app <- ExploreModelMatrix(sampleData = coldata, 
                         designFormula = ~0+ vector + LPS + vector:LPS)
#if (interactive()) shiny::runApp(app)

#See script "checking_ways_of_specifying_contrasts" for more details
#to understand how the design and contrasts work

#differential expression analysis

setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/DESeq2_results")

sink(paste0("DESeq2.log"))
dds <- DESeq(ddsMat)
resultsNames(dds)

comparison_list <- list(p30vsEV=c("vector","p30","EV"),
                        p42vsEV=c("vector","p30","p42"),
                        p30vsp42=c("vector","p30","p42"),
                        p30LPSvsEVLPS=c(-1,1,0,0,1,0),
                        p42LPSvsEVLPS=c(-1,0,1,0,0,1),
                        p30LPSvsp42LPS=c(0,1,-1,0,1,-1),
                        EVLPSvsEV =c(0,0,0,1,0,0),
                        p30LPSvsp30=c(0,0,0,1,1,0),
                        p42LPSvsp42=c(0,0,0,1,0,1))


for(comp in names(comparison_list)){
  

  Group <- t(metadata%>%filter(sample==gene)%>%select(-1))
  Group[Group==0] <- paste0(gene,"_WT")
  Group[Group==1] <- paste0(gene,"_MU")
  Group <- as.factor(Group)
  dds <- DESeqDataSetFromMatrix(round(rsem_counts), DataFrame(Group), ~ 0+ Group)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds$Group <- relevel(dds$Group, ref =paste0(gene,"_WT"))
  dds <- DESeq(dds)
  res <- results(dds,alpha = 0.05)
  res <- res[order(res$padj),]
  summary(res)
  
  write.csv(res,paste0("results/DESeq2/DESeq2_results_",gene,".csv"),quote = F)
  rank <- as.data.frame(res[order(res$stat,decreasing = T),])
  rank$gene <- rownames(rank)
  write.table(rank[,c("gene","stat")],paste0("results/DESeq2/rank_files/",gene,"_DESeq2.rnk"),quote = F,col.names = F,row.names = F,sep = "\t")
  sink()
}

results(dds, contrast = c("vector","p30","EV"))
head(p30_vs_EV)
summary(p30_vs_EV)

#P30 LPS vs p42 LPS


summary(p30LPS_vs_p42LPS)
head(p30LPS_vs_p42LPS[order(p30LPS_vs_p42LPS$padj),])


res <- results(dds,alpha = 0.05)
res <- res[order(res$padj),]
summary(res)

write.csv(res,paste0("results/DESeq2/DESeq2_results_",gene,".csv"),quote = F)

res <- results(dds)
head(res)
summary(res)
resultsNames(dds)

#desired comparisons:

#P30 vs EV
p30_vs_EV <- results(dds, contrast = c("vector","p30","EV"))
head(p30_vs_EV)
summary(p30_vs_EV)

#P42 vs EV
p42_vs_EV <- results(dds, contrast = c("vector","p42","EV"))
head(p42_vs_EV)
summary(p42_vs_EV)

# P30 vs p42 
p30vsp42 <- results(dds, contrast=c("vector","p30","p42"))
p30vsp42
summary(p30vsp42)

#p42_LPS vs EV_LPS 
p42LPS_vs_EVLPS <- results(dds,c(-1,0,1,0,0,1))
summary(p42LPS_vs_EVLPS)
#p42LPS_vs_EVLPS_alt <- results(dds,contrast = list(c("vectorp42" ,"vectorp42.LPSTRUE"),"vectorEV"))
#summary(p42LPS_vs_EVLPS_alt)


# P30_LPS vs EV_LPS 
p30LPS_vs_EVLPS <- results(dds,c(-1,1,0,0,1,0))
summary(p30LPS_vs_EVLPS)
#p30LPS_vs_EVLPS_alt <- results(dds,contrast = list(c("vectorp30" ,"vectorp30.LPSTRUE"),"vectorEV"))
#summary(p30LPS_vs_EVLPS_alt)


# EV_LPS vs EV 
EVLPS_vs_EV <- results(dds,c(0,0,0,1,0,0))
summary(EVLPS_vs_EV)

# P30_LPS vs P30 
p30LPS_vs_p30 <- results(dds,c(0,0,0,1,1,0))
summary(p30LPS_vs_p30)

# P42_LPS vs p42 
p42LPS_vs_p42 <- results(dds,c(0,0,0,1,0,1))
summary(p42LPS_vs_p42)

#P30 LPS vs p42 LPS
p30LPS_vs_p42LPS=results(dds,contrast = c(0,1,-1,0,1,-1))

summary(p30LPS_vs_p42LPS)
head(p30LPS_vs_p42LPS[order(p30LPS_vs_p42LPS$padj),])


#Now make sense of the gene names
library(rtracklayer)
geneannot <- readGFF("../../../references/mouse/gencode/gencode.vM27.annotation.gtf")
geneannot_genes <- geneannot%>%filter(type=="gene")
length(unique(geneannot_genes$gene_id))
length(unique(geneannot_genes$gene_name))
length(unique(paste(geneannot_genes$gene_id,geneannot_genes$gene_name)))

duplicated_gene_names <- unique(geneannot_genes$gene_name[duplicated(geneannot_genes$gene_name)])
length(duplicated_gene_names)
duplicated_gene_names
geneannot_genes_dup <- geneannot_genes%>%filter(gene_name%in%duplicated_gene_names)
length(unique(paste(geneannot_genes_dup$gene_type,geneannot_genes_dup$gene_name)))
table(geneannot_genes_dup$gene_type)

test=read_rds("~/Downloads/Mm.c1.all.v7.1.entrez.rds")
test2=read_rds("~/Downloads/Mm.c5.all.v7.1.entrez.rds")
Hallmarks_mouse <- read_rds("~/Downloads/Mm.h.all.v7.1.entrez.rds")
mm_ncbi_geneinfo=fread("~/Downloads/Mus_musculus.gene_info")

table(geneannot_genes$gene_name%in%mm_ncbi_geneinfo$Symbol)
table(unique(unlist(Hallmarks_mouse))%in%mm_ncbi_geneinfo$GeneID)

allsynonyms <- unlist(sapply(mm_ncbi_geneinfo$Synonyms,function(x)strsplit(x,split = "\\|")))

Hallmarks_mouse_unique_ids <- unique(unlist(Hallmarks_mouse))
Hallmarks_mouse_unique_Symbols <- mm_ncbi_geneinfo$Symbol[match(Hallmarks_mouse_unique_ids,
                                                                mm_ncbi_geneinfo$GeneID,nomatch = 0)]
length(unique(Hallmarks_mouse_unique_Symbols))
table(Hallmarks_mouse_unique_Symbols%in%geneannot_genes$gene_name)


table(geneannot_genes$gene_name%in%c(allsynonyms,mm_ncbi_geneinfo$Symbol,mm_ncbi_geneinfo$Symbol_from_nomenclature_authority))
geneannot_genes$gene_name[!geneannot_genes$gene_name%in%c(allsynonyms,mm_ncbi_geneinfo$Symbol,mm_ncbi_geneinfo$Symbol_from_nomenclature_authority)]

#convert ids in Hallmark genes to Symbols or vice versa?
head(DESeq2::counts(dds))

conversion_file <- fread("~/share/Cuartero Group/CUARTERO GROUP/references/mouse/Ensembl_to_NCBI_id_conversion.tsv",data.table = F)
table(geneannot_genes$gene_id%in%conversion_file$`Gene stable ID version`)
tmplog <- !duplicated(conversion_file[,c("Gene stable ID version","NCBI gene (formerly Entrezgene) ID")])
conversion_file_undup <- conversion_file[tmplog,]
conversion_file_undup_ncbiIDavail <- conversion_file_undup[!is.na(conversion_file_undup$`NCBI gene (formerly Entrezgene) ID`),]
length(unique(conversion_file_undup_ncbiIDavail$`Gene stable ID version`))
length(unique(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`))

duplicatedtmp <- unique(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`[duplicated(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`)])
table(geneannot_genes$gene_id%in%conversion_file_undup_ncbiIDavail$`Gene stable ID version`)
table(conversion_file_undup_ncbiIDavail$`Gene stable ID version`%in%geneannot_genes$gene_id)

#note that some Ensembl ids match with more than one NCBI id
#then try to select the one that is on the gene set db
table(Hallmarks_mouse_unique_ids%in%conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`)
table(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`%in%Hallmarks_mouse_unique_ids)
#sort table in a way that those ids that are in Hallmarks appear first:
logivec <- conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`%in%Hallmarks_mouse_unique_ids
head(order(logivec,decreasing = T))
head(logivec,14)

conversion_file_undup_ncbiIDavail_sorted <- 
  conversion_file_undup_ncbiIDavail[order(logivec,decreasing = T),]

conversion_file_undup_ncbiIDavail_sorted$`NCBI gene (formerly Entrezgene) ID`%in%Hallmarks_mouse_unique_ids

geneannot_genes$ncbi_GeneID <- conversion_file_undup_ncbiIDavail_sorted$`NCBI gene (formerly Entrezgene) ID`[match(geneannot_genes$gene_id,
                                                                                                                   conversion_file_undup_ncbiIDavail_sorted$`Gene stable ID version`)]
rownames.to.change <- rownames(countdata)[rownames(countdata)%in%geneannot_genes$gene_id[!is.na(geneannot_genes$ncbi_GeneID)]] 
names.to.change.to <- geneannot_genes$ncbi_GeneID[match(rownames.to.change,geneannot_genes$gene_id,nomatch = 0)]

length(unique(names.to.change.to))
names.to.change.to.dup <- unique(names.to.change.to[duplicated(names.to.change.to)])
table(names.to.change.to.dup%in%Hallmarks_mouse_unique_ids)
table(Hallmarks_mouse_unique_ids%in%names.to.change.to)
Hallmarks_mouse_unique_ids_not.in.gencode <- Hallmarks_mouse_unique_ids[!Hallmarks_mouse_unique_ids%in%geneannot_genes$ncbi_GeneID]
View(conversion_file[conversion_file$`NCBI gene (formerly Entrezgene) ID`%in%Hallmarks_mouse_unique_ids_not.in.gencode,])

#some genes are missing lost in translation, and some
#are duplicated...I will just remove those

rownames.to.change <- rownames.to.change[!duplicated(names.to.change.to)]
names.to.change.to <- names.to.change.to[!duplicated(names.to.change.to)]

rownames(countdata)[rownames(countdata)%in%rownames.to.change] <- names.to.change.to

table(rownames(countdata)%in%Hallmarks_mouse_unique_ids)
table(Hallmarks_mouse_unique_ids%in%rownames(countdata))

#there are 128 Hallmark genes that are missing in the countdata

#dseq2 tables
#fpkm gene_name