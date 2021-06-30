#Now make sense of the gene names
library(rtracklayer)
setwd("share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/scripts/")
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

c2=read_rds("~/Downloads/Mm.c2.all.v7.1.entrez.rds")
c3=read_rds("~/Downloads/Mm.c3.all.v7.1.entrez.rds")
c7=read_rds("~/Downloads/Mm.c7.all.v7.1.entrez.rds")
#Hallmarks_mouse <- read_rds("~/Downloads/Mm.h.all.v7.1.entrez.rds")
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