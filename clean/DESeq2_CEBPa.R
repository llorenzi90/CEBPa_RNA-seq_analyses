library(DESeq2)
setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/tsv_files/")
#setwd("/Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results")
sink(paste0("DESeq2.log"))

countdata=read.table("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.txt",header = T)

#read conversion table for rank files (mouse gene sets have NCBI ids)
conversion_ENSMBL_NCBI <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/references/mouse/updated_conversion_table_ENSEMBL_NCBI_ids.csv")
conversion_ENSMBL_NCBI <- conversion_ENSMBL_NCBI[!is.na(conversion_ENSMBL_NCBI$ncbi_GeneID),]

rownames(countdata) <- countdata$gene_name
countdata <- countdata[,-1]
name=colnames(countdata)
tr=strsplit(name,split = "\\.")
vector=sapply(tr, function(x)return(x[1]))
rep=sapply(tr, function(x)return(x[length(x)]))
LPS=sapply(tr, function(x)return(any(grepl("LPS",x))))
coldata=data.frame(name,vector,LPS,rep)
rownames(coldata) <- coldata$name
#table(rownames(coldata)==colnames(countdata))
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

keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]
nrow(ddsMat)

ddsMat <- estimateSizeFactors(ddsMat)

dds <- DESeq(ddsMat)
resultsNames(dds)

comparison_list <- list(p30vsEV=c("vector","p30","EV"),
                        p42vsEV=c("vector","p42","EV"),
                        p30vsp42=c("vector","p30","p42"),
                        p30LPSvsEVLPS=c(-1,1,0,0,1,0),
                        p42LPSvsEVLPS=c(-1,0,1,0,0,1),
                        p30LPSvsp42LPS=c(0,1,-1,0,1,-1),
                        EVLPSvsEV =c(0,0,0,1,0,0),
                        p30LPSvsp30=c(0,0,0,1,1,0),
                        p42LPSvsp42=c(0,0,0,1,0,1))

results_list <- list()

for(comp in names(comparison_list)){
  cat("\n")
  print(comp)
  res <- results(dds, contrast = comparison_list[[comp]],alpha = 0.05)
  summary(res)
  res <- res[order(res$padj),]
  results_list[[comp]] <- res
  
  write.table(res,paste0(comp,"_DESeq_results.txt"),quote = F,sep = "\t")
  rank <- as.data.frame(res[order(res$stat,decreasing = T),])
  rank$gene <- rownames(rank)
  rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$gene_name] <- 
    conversion_ENSMBL_NCBI$ncbi_GeneID[match(rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$gene_name],
                                             conversion_ENSMBL_NCBI$gene_name)]
  write.table(rank[,c("gene","stat")],paste0("rank_files/",comp,"_DESeq2.rnk"),quote = F,col.names = F,row.names = F,sep = "\t")

}

sink()
