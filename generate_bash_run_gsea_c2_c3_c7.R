## ---------------------------
##
##
## Purpose of script: create script to run GSEA for all samples with DESeq2 rank files
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-05-29
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

#command example: gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.3.symbols.gmt -collapse Remap_Only -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/llorenzi/Nextcloud/TCGA/LL_analyses/TCGA_data_analysis/results/limma_voom/rank_files/ASXL1_limma_voom.rnk -scoring_scheme weighted -rpt_label my_analysis -chip /Users/llorenzi/Nextcloud/TCGA/LL_analyses/TCGA_data_analysis/data/GSEA/match_gene_sets.chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/llorenzi/gsea_home/output/mar31
#cmm1 <- "./gsea-cli.sh GSEAPreranked -gmx /Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/Mm.h.all.v7.1.entrez.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/rank_files/"

gene_sets <- c("Mm.c2.all.v7.1.entrez.gmt",
               "Mm.c3.all.v7.1.entrez.gmt",
               "Mm.c7.all.v7.1.entrez.gmt")
setwd("/home")

cmm1 <- './gsea-cli.sh GSEAPreranked -gmx ../share/Cuartero\\ Group/CUARTERO\\ GROUP/references/mouse/'
cmm1.1 <- ' -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ../share/Cuartero\\ Group/CUARTERO\\ GROUP/CBPa/RNA-seq/DESeq2_results/rank_files/'
cmm2 <- " -scoring_scheme weighted -rpt_label "
#cmm3 <- " -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/llorenzi/Nextcloud/CBPa_RNA_seq/GESEA_results"
cmm3 <- " -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ../share/Cuartero\\ Group/CUARTERO\\ GROUP/CBPa/RNA-seq/GSEA_results/"

commands <- c()
for (geneset in gene_sets) {
  dir.create(paste0("~/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/GSEA_results/",geneset))
  #for(f in list.files("/Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/rank_files/")){
  for(f in list.files("~/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/DESeq2_results/rank_files/")){
    
    nam <- gsub("_DESeq2.rnk","",f)
    commands <- c(commands, paste0(cmm1,geneset,cmm1.1,f,cmm2,nam,cmm3,geneset))
  }
  
  #write(commands,"/Users/llorenzi/Nextcloud/CBPa_RNA_seq/run_GseaPreranked_CBPa_DESeq2.sh")
  write(commands,"~/share/Cuartero Group/CUARTERO GROUP/CBPa/RNA-seq/run_GseaPreranked_CBPa_DESeq2_new_gene_sets.sh")
  
  #run this script from /Users/llorenzi/GSEA_4.1.0 folder
    
}
