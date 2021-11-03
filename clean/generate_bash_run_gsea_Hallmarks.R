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

#GESEA version 4.1.0
#~/software/GSEA_Linux_4.1.0/
cmm1 <- './gsea-cli.sh GSEAPreranked -gmx ../../share/Cuartero\\ Group/CUARTERO\\ GROUP/references/mouse/Mm.h.all.v7.1.entrez.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ../../share/Cuartero\\ Group/CUARTERO\\ GROUP/CEBPa/RNA-seq/results/DESeq2_results/rank_files/'
cmm2 <- " -scoring_scheme weighted -rpt_label "
cmm3 <- " -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ../../share/Cuartero\\ Group/CUARTERO\\ GROUP/CEBPa/RNA-seq/results/GSEA_results/Hallmarks"

commands <- c()
for(f in list.files("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/results/DESeq2_results/rank_files/")){
  nam <- gsub("_DESeq2.rnk","",f)
  commands <- c(commands, paste0(cmm1,f,cmm2,nam,cmm3))
}

write(commands,"~/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/scripts/clean/run_GseaPreranked_CEBPa_DESeq2.sh")

#run this script from ~/software/GSEA_Linux_4.1.0/