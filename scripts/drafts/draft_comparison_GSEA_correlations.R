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
  
  
    

sessionInfo()
  
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.4.0      data.table_1.14.0 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.6      
# [6] purrr_0.3.4       readr_1.4.0       tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.4    
# [11] tidyverse_1.3.1  
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.1 splines_4.0.3    haven_2.4.1      lattice_0.20-44  carData_3.0-4   
# [6] colorspace_2.0-1 vctrs_0.3.8      generics_0.1.0   mgcv_1.8-36      utf8_1.2.1      
# [11] rlang_0.4.11     pillar_1.6.1     foreign_0.8-81   glue_1.4.2       withr_2.4.2     
# [16] DBI_1.1.1        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1     lifecycle_1.0.0 
# [21] ggsignif_0.6.2   munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 zip_2.2.0       
# [26] rvest_1.0.0      rio_0.5.26       labeling_0.4.2   curl_4.3.1       fansi_0.5.0     
# [31] broom_0.7.7      Rcpp_1.0.6       polynom_1.4-0    scales_1.1.1     backports_1.2.1 
# [36] jsonlite_1.7.2   abind_1.4-5      farver_2.1.0     fs_1.5.0         hms_1.1.0       
# [41] digest_0.6.27    openxlsx_4.2.4   stringi_1.6.2    rstatix_0.7.0    grid_4.0.3      
# [46] cli_2.5.0        tools_4.0.3      magrittr_2.0.1   car_3.0-10       crayon_1.4.1    
# [51] pkgconfig_2.0.3  ellipsis_0.3.2   Matrix_1.3-4     xml2_1.3.2       reprex_2.0.0    
# [56] lubridate_1.7.10 assertthat_0.2.1 httr_1.4.2       rstudioapi_0.13  R6_2.5.0        
# [61] nlme_3.1-152     compiler_4.0.3  
#     
# 
# 
# 
