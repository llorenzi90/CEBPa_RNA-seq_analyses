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
library(ggpubr)
setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/GSEA_results/Comparisons_TCGA")
setwd("~/CEBPa/results/GSEA_results/Comparisons_TCGA/")

pltf <- function(da,tit="",colby="significant",showcor=T,meth="pearson"){
  g <- ggplot(da,aes(x=`HPC-7`,y=TCGA)) +
    geom_point(aes_string(color=colby))+
    geom_smooth(method = lm)+
    ggtitle(paste0(geneset, " NES correlation\n(",tit,")"))+ 
    theme_classic()
  if(showcor){
    print(g+stat_cor(method = meth))
  }else print(g+ stat_regline_equation(label.y = 0.7*max(g$data$TCGA,na.rm = T), aes(label = ..rr.label..)))
}


####Preset parameters:

##RUN 1
wd <- "~/CEBPa/results/GSEA_results/Comparisons_TCGA/"
setwd(wd)
dir.create("correlation_plots")  
wd <- paste0(wd,"correlation_plots")
wd
#showcor=T
showcor=F
meth="pearson"


###RUN 2
wd <- "~/CEBPa/results/GSEA_results/Comparisons_TCGA/"
setwd(wd)
dir.create("correlation_plots_pearson")  
wd <- paste0(wd,"correlation_plots_pearson")
wd
showcor=T
#showcor=F
meth="pearson"

###RUN 3
wd <- "~/CEBPa/results/GSEA_results/Comparisons_TCGA/"
setwd(wd)
dir.create("correlation_plots_spearman")  
wd <- paste0(wd,"correlation_plots_spearman")
wd
showcor=T
#showcor=F
meth="spearman"

for (geneset in c("Hallmarks","c2","c3","c7")) {
  comparison_mouse <- "p30vsp42"
  comparison_TCGA <- "CEBPA"
  #datadir_mouse <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/GSEA_results"
  datadir_mouse <- "~/CEBPa/results/GSEA_results"
  
  #datadir_TCGA <- "/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/TCGA/LL_analyses/results/DESeq2_GSEA"
  datadir_TCGA <- "/Users/llorenzi/Documents/TCGA/LL_analyses/results/DESeq2_GSEA"
  
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
  
  plot(data_all$NES[data_all$X12=="HPC-7"],
       data_all$NES[data_all$X12=="TCGA"])
  
  colnames(data_all)[12] <- "dataset"
  
  
  for(qval in c(0.05,0.01)){
    setwd(wd)
    outdir <- paste0("qval_",qval)
    dir.create(outdir)
    setwd(outdir)
    
    #reorder data
    df_qval <- spread(data_all[,c("NAME","FDR q-val","dataset")],key=dataset,value = `FDR q-val`)
    
    df_qval <- df_qval %>% mutate(significant = case_when(`HPC-7`<qval & TCGA>=qval ~ "HPC-7-only",
                                                           `HPC-7`>=qval & TCGA<qval ~ "TCGA-only",
                                                           `HPC-7`<qval & TCGA<qval ~ "both",
                                                           `HPC-7`>=qval & TCGA>=qval ~ "none")) 
    
    df_corr_mut <- spread(data = data_all[,c("NAME","dataset","NES")],key = dataset,value = NES)
    df_corr_mut <- as.data.frame(df_corr_mut)
    df_corr_mut$significant <- df_qval$significant[match(df_corr_mut$NAME,
                                                           df_qval$NAME)]
    
    #make plot with all datapoints (i.e not only those that are
    #significant in both datasets)
    ####save plots to a single file####
    
    pdf(paste0(geneset,"_correlation_plots.pdf"),onefile = T)
    ##plot 1 All genesets, colour by significance
    da <- df_corr_mut  
    tit <- "all genesets"
    pltf(da,tit,showcor = showcor,meth = meth)
    
    ##plot 2 At least one significant, colour by significance
    da <- df_corr_mut%>%filter(significant!="none")
    tit <- "significant in at least one"
    pltf(da,tit,showcor = showcor,meth = meth)
    
    #plot 3 significant in both, colour by significance
    
    da <- df_corr_mut%>%filter(significant=="both")
    tit <- "significant in both"
    pltf(da,tit,showcor = showcor,meth = meth)
    
    #> `geom_smooth()` using formula 'y ~ x'
    #plot 4 significant in both, colour by geneset
    da <- df_corr_mut%>%filter(significant=="both")
    tit <- "significant in both"
    
    #change names to make them shorter
    da$NAME <-  abbreviate(gsub("_",".",gsub("HALLMARK_","",da$NAME)),minlength = 18)
    
    da$NAME <- factor(da$NAME,   levels = unique(da$NAME))
    
    pltf(da,tit,colby = "NAME",showcor = showcor,meth = meth)
    dev.off()
    
    
    }
  }
  
  
    

  
  
  
    



