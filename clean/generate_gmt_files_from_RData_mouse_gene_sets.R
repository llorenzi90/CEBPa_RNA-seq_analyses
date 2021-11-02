##read hallmark RData files gene sets for mouse:
#source of GSEA mouse hallmarks: http://bioinf.wehi.edu.au/MSigDB/

# Background
# The Molecular Signatures Database (MSigDB) is an important resource created and maintained by the Broad Institute. The gene sets contained in the MSigDB are from a wide variety of sources and consist of human genes identified either by NCBI GeneID or by gene symbol. Our work at the WEHI predominately uses mouse models of human disease. To facilitate use of the MSigDB in our work, we have created a pure mouse version of the MSigDB by mapping all sets to mouse orthologs. A pure human version is also provided.
# 
# Procedure
# For human, gmt files downloaded from the MSigDB were converted to R lists and saved in RDS format.
# 
# The mouse C1 positional gene set collection was created from the NCBI gene information file Mus_musculus.gene_info.gz downloaded from the NCBI ftp site. Cytobands were identified from the map_location column. The positional collection provides GeneIDs for the genes in each cytoband.
# 
# The C5 gene ontology collections were created from the Bioconductor organism package org.Mm.eg.db using the GeneID to GO Term mappings provided by the egGO2ALLEGS Bimap. This ensures that the GO Term hierarchy is respected: any GeneID associated with a child (more specific) GO Term is also included in any parent (more general) GO Term that is an ancestor for the original Term.
# 
# All other mouse collections were created by mapping the corresponding human collection to mouse orthologs using HGNC Comparison of Orthology Predictions (HCOP). The HCOP tool integrates orthology assertions predicted by eggNOG, Ensembl Compara, HGNC, HomoloGene, Inparanoid, NCBI Gene Orthology, OMA, OrthoDB, OrthoMCL, Panther, PhylomeDB, PomBase, TreeFam and ZFIN. It includes non-coding as well a protein-coding genes.
# 
# Current Version
# RDS files for MSigDB Version 7.1, based on MSigDB version 7.1, org.Mm.eg.db version 3.11.4 and NCBI gene information file dated 12 June 2020.

Hallmarks_mouse <- read_rds("~/Downloads/Mm.h.all.v7.1.entrez.rds")
c2_mouse <- read_rds("~/Downloads/Mm.c2.all.v7.1.entrez.rds")
c3_mouse <- read_rds("~/Downloads/Mm.c3.all.v7.1.entrez.rds")
c7_mouse <- read_rds("~/Downloads/Mm.c7.all.v7.1.entrez.rds")
gene_sets_mouse_unique_ids <- unique(unlist(sapply(list(Hallmarks_mouse,
                                                        c2_mouse,
                                                        c3_mouse,
                                                        c7_mouse),function(x)unique(unlist(x)))))
head(gene_sets_mouse_unique_ids)

##### generate gmt files ######
lapply(names(Hallmarks_mouse),function(x){
  write(paste(x,NA,paste(Hallmarks_mouse[[x]],collapse = "\t"),sep = "\t"),"Mm.h.all.v7.1.entrez.gmt",append = T)
})

new_gene_sets_list <- list(c2_mouse,c3_mouse,c7_mouse)
names(new_gene_sets_list) <- c("Mm.c2.all.v7.1.entrez",
                               "Mm.c3.all.v7.1.entrez",
                               "Mm.c7.all.v7.1.entrez")

for (nam in names(new_gene_sets_list)) {
  tmp <- new_gene_sets_list[[nam]]
  lapply(names(tmp),function(x){
    write(paste(x,NA,paste(tmp[[x]],collapse = "\t"),sep = "\t"),paste0("~/share/Cuartero Group/CUARTERO GROUP/references/mouse/",nam,".gmt"),append = T)
  })
}
