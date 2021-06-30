## ---------------------------
##
##
## Purpose of script: plot gene set enrichment analysis 
##                    to compare HPC-7 mouse cells 
##                    (which do not endogenously express 
##                    C/EBPa) transfected with either the 
##                    p42 isoform or the p30 isoform 
##                    vs TCGA patients with or without
##                    CBPa mutations
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-06-14
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
setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/GSEA_results/Comparisons_TCGA")
## ---------------------------
#geneset <- "Hallmarks"
#geneset <- "c2"
#geneset <- "c3"
#geneset <- "c7"
# custom_name_abbreviation <- function(str,
#                                      flength=15,
#                                      split="_",
#                                      minlen=3,
#                                      maxterms=4){
#   x <- strsplit(str,split = split)[[1]]
#   
#   tmp=1
#   while(length(x)>4){
#     if(!all(nchar(x)<minlen)& any(nchar(x)<minlen))x <- x[-which(nchar(x)<minlen)[1]]
#     else  if((tmp %% 2) != 0) {
#       x <- x[-length(x)]
#       tmp <- tmp +1}else x <- x[-1]
#       
#   }
#   l <- length(x)
#   minlelement <- floor(flength/l)
#   lengths <- rep(minlelement,l)
#   resto <- flength - minlelement*l
#   if(resto!=0)  lengths[1:resto] <- lengths[1:resto] +1
#   
#   paste(sapply(seq_along(x), 
#                function(z)abbreviate(x[z],strict = T,minlength = lengths[z])),collapse = ".")
# }
#test <- sapply(as.character(data_all$NAME),custom_name_abbreviation)

for (geneset in c("Hallmarks","c2","c3","c7")) {
  comparison_mouse <- "p30vsp42"
  comparison_TCGA <- "CEBPA"
  datadir_mouse <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/GSEA_results"
  datadir_TCGA <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/TCGA/LL_analyses/results/DESeq2_GSEA"
  
  datadir_mouse <- list.files(datadir_mouse,full.names = T,pattern = geneset )
  datadir_mouse <- list.files(datadir_mouse,full.names = T,pattern = comparison_mouse)
  datapaths_mouse <- list.files(datadir_mouse,full.names = T,pattern = "gsea_report_for.*tsv")
  
  datadir_TCGA <- list.files(datadir_TCGA,full.names = T,pattern = geneset )
  datadir_TCGA <- list.files(datadir_TCGA,full.names = T,pattern = comparison_TCGA)
  datapaths_TCGA <- list.files(datadir_TCGA,full.names = T,pattern = "gsea_report_for.*tsv")
  
  data_mouse <- map_dfr(datapaths_mouse,function(x)read_tsv(x)%>%mutate(NES=as.numeric(NES),
                                                                        `NOM p-val`=as.numeric(`NOM p-val`)))
  
  data_TCGA <- map_dfr(datapaths_TCGA,function(x)read_tsv(x)%>%mutate(NES=as.numeric(NES),
                                                                      `NOM p-val`=as.numeric(`NOM p-val`)))
  
  table(data_mouse$NAME%in%data_TCGA$NAME)
  table(data_TCGA$NAME%in%data_mouse$NAME)
  
  data_TCGA <- data_TCGA %>% slice(match(data_mouse$NAME,data_TCGA$NAME,nomatch = 0))
  data_mouse <- data_mouse %>% filter(data_mouse$NAME%in%data_TCGA$NAME)
  
  table(data_mouse$NAME==data_TCGA$NAME)
  
  #plot(data_mouse$NES,data_TCGA$NES)
  #cor(data_mouse$NES,data_TCGA$NES,method = "spearman")
  
  
  data_all <- bind_rows(data_mouse %>% mutate(X12="HPC-7"),
                        data_TCGA %>% mutate(X12="TCGA"))
  
  
  colnames(data_all)[12] <- "dataset"
  for(qval in c(0.05,0.01)){
    for (sortmeth in c("qval","NES")) {
      if(sortmeth=="qval"){
        data_all_filtered <- data_all%>%group_by(NAME) %>%
          filter(all(`FDR q-val`<qval)) %>%
          mutate(mi=min(`FDR q-val`))%>%
          mutate(mx=max(`FDR q-val`))%>%
          arrange(mi,mx,NAME)%>%
          select(-mx,-mi)}else {data_all_filtered <- data_all%>%group_by(NAME) %>%
              filter(all(`FDR q-val`<qval)) %>%
              mutate(mi=min(NES))%>%
              mutate(mx=max(NES))%>%
              arrange(mi,mx,NAME)%>%
              select(-mx,-mi)
          }
      data_all_filtered$NAME <- abbreviate(gsub("_",".",gsub("HALLMARK_","",data_all_filtered$NAME)),minlength = 18)
      
        data_all_filtered$NAME <- factor(data_all_filtered$NAME,
                                         levels = unique(data_all_filtered$NAME))
        g <- ggplot(data_all_filtered,aes(x=NES,y=NAME)) +
          geom_point(aes(color=dataset,size=`FDR q-val`))+
          scale_size(trans = 'reverse')+ 
          ggtitle(comparison_mouse)+ 
          ylab(geneset)+
          geom_vline(xintercept = 0, linetype="dashed")+
          theme_classic()
        
        pdf(paste0("test_",geneset,"_filtered_min_qval_",qval,"_sorted_by_",sortmeth,".pdf"))
        print(g)
        dev.off()
      }
    }
   
    
    
    
  }
  
  
  
    



