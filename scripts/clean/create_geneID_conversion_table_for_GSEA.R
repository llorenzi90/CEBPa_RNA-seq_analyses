library(rtracklayer)
geneannot <- readGFF("~/share/Cuartero Group/CUARTERO GROUP/references/mouse/gencode/gencode.vM27.annotation.gtf")
geneannot_genes <- geneannot%>%filter(type=="gene")
length(unique(geneannot_genes$gene_id))
length(unique(geneannot_genes$gene_name))
length(unique(paste(geneannot_genes$gene_id,geneannot_genes$gene_name)))
rownames_gene_id_counts <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/htseq_count/data/all_samples_gencodevM27.counts.csv")[,1]
rownames_gene_name_counts <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv")[,1]

table(geneannot_genes$gene_name[match(rownames_gene_id_counts,geneannot_genes$gene_id)]==rownames_gene_name_counts)
rownames_gene_name_counts[!geneannot_genes$gene_name[match(rownames_gene_id_counts,geneannot_genes$gene_id)]==rownames_gene_name_counts]

##read hallmark gene sets for mouse:
#source of GSEA mouse hallmarks: http://bioinf.wehi.edu.au/MSigDB/
Hallmarks_mouse <- read_rds("~/Downloads/Mm.h.all.v7.1.entrez.rds")
c2_mouse <- read_rds("~/Downloads/Mm.c2.all.v7.1.entrez.rds")
c3_mouse <- read_rds("~/Downloads/Mm.c3.all.v7.1.entrez.rds")
c7_mouse <- read_rds("~/Downloads/Mm.c7.all.v7.1.entrez.rds")
gene_sets_mouse_unique_ids <- unique(unlist(sapply(list(Hallmarks_mouse,
                                          c2_mouse,
                                          c3_mouse,
                                          c7_mouse),function(x)unique(unlist(x)))))
head(gene_sets_mouse_unique_ids)
#these are NCBI ids

#read conversion table generetad with biomart:
#conversion_file <- fread("~/Downloads/mouse/Ensembl_to_NCBI_id_conversion.tsv",data.table = F)
conversion_file <- fread("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/references/mouse/Ensembl_to_NCBI_id_conversion.tsv",data.table = F)
#all ENSEMBL ids in our countdata are present in this table
table(geneannot_genes$gene_id%in%conversion_file$`Gene stable ID version`)
table(gene_sets_mouse_unique_ids%in%conversion_file$`NCBI gene (formerly Entrezgene) ID`)
20462/length(gene_sets_mouse_unique_ids)*100
#95% of the NCBI ids for hallmarks are present

#remove redundant rows 
tmplogic <- !duplicated(conversion_file[,c("Gene stable ID version","NCBI gene (formerly Entrezgene) ID")])
table(tmplogic)
conversion_file_undup <- conversion_file[tmplogic,]
#remove rows with no ncbi id (NA)
conversion_file_undup_ncbiIDavail <- conversion_file_undup[!is.na(conversion_file_undup$`NCBI gene (formerly Entrezgene) ID`),]

nrow(conversion_file_undup_ncbiIDavail)
length(unique(conversion_file_undup_ncbiIDavail$`Gene stable ID version`))
length(unique(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`))
#there are both ENSEMBL and NCBI ids duplicated 
#I want to remove duplicated ids 
#prioritizing leaving those that are present in hallmarks
#duplicated ncbi ids:
duplicatedtmp <- unique(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`[duplicated(conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`)])
table(geneannot_genes$gene_id%in%conversion_file_undup_ncbiIDavail$`Gene stable ID version`)
table(conversion_file_undup_ncbiIDavail$`Gene stable ID version`%in%geneannot_genes$gene_id)

#note that some Ensembl ids match more than one NCBI id
#then try to select one that is on the gene set db
#sort table in a way that those ids that are in Hallmarks appear first:
logivec <- conversion_file_undup_ncbiIDavail$`NCBI gene (formerly Entrezgene) ID`%in%gene_sets_mouse_unique_ids
head(order(logivec,decreasing = T))
head(logivec,14)

conversion_file_undup_ncbiIDavail_sorted <- 
  conversion_file_undup_ncbiIDavail[order(logivec,decreasing = T),]

conversion_file_undup_ncbiIDavail_sorted$`NCBI gene (formerly Entrezgene) ID`%in%gene_sets_mouse_unique_ids

geneannot_genes$ncbi_GeneID <- conversion_file_undup_ncbiIDavail_sorted$`NCBI gene (formerly Entrezgene) ID`[match(geneannot_genes$gene_id,
                                                                                                                   conversion_file_undup_ncbiIDavail_sorted$`Gene stable ID version`)]
length(unique(geneannot_genes$ncbi_GeneID))
table(gene_sets_mouse_unique_ids%in%geneannot_genes$ncbi_GeneID)
table(!is.na(geneannot_genes$ncbi_GeneID))
table(duplicated(geneannot_genes$ncbi_GeneID[!is.na(geneannot_genes$ncbi_GeneID)]))

length(unique(geneannot_genes$ncbi_GeneID[duplicated(geneannot_genes$ncbi_GeneID)&!is.na(geneannot_genes$ncbi_GeneID)]))
#note that 47 ncbi ids are duplicated 
duplicated_ncbi_ids <- unique(geneannot_genes$ncbi_GeneID[duplicated(geneannot_genes$ncbi_GeneID)&!is.na(geneannot_genes$ncbi_GeneID)])
table(duplicated_ncbi_ids%in%gene_sets_mouse_unique_ids)
#there are 34 that are also in hallmarks
table(gene_sets_mouse_unique_ids%in%duplicated_ncbi_ids|(!gene_sets_mouse_unique_ids%in%geneannot_genes$ncbi_GeneID))
#in total there are 1139 genes that are either
#ambiguous relative to ENSEMBL ids (match multiple ids) 
#or have no match at all
#this represents 
1139/length(gene_sets_mouse_unique_ids)*100
#5.3% of all genes, therefore I will just remove those
#make a conversion table with those genes
#that do have a match in Ensembl names/ids:
rownames <- cbind(rownames_gene_id_counts,rownames_gene_name_counts)
colnames(rownames) <- c("gene_id","gene_name")
rownames <- as.data.frame(rownames)

#selected ncbi ids:
sel_geneannot_genes <- geneannot_genes[!is.na(geneannot_genes$ncbi_GeneID)&!duplicated(geneannot_genes$ncbi_GeneID),]
table(gene_sets_mouse_unique_ids%in%sel_geneannot_genes$ncbi_GeneID)

rownames$ncbi_GeneID <- sel_geneannot_genes$ncbi_GeneID[match(rownames$gene_id,sel_geneannot_genes$gene_id)]
write.csv(rownames,"updated_conversion_table_ENSEMBL_NCBI_ids.csv",row.names = F,quote = F)

