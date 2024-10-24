## The role of CEBPa in the inflammatory response 
HPC-7 cells (which do not endogenously express C/EBPa) were transfected with plasmids containing either the p42 isoform or the p30 isoform. P42 and p30 are bonded to ERT2 (estrogen receptor), so that it localises in the cytoplasm until its induction with tamoxifen (ligand), which will translocate to the nucleus and only then p42 and p30 will be able to function. This fusion was performed to avoid undesired differentiative effects derived from the presence of the p42 isoform in the nucleus. We have also observed that AML patients carrying this C/EBPa mutation have a reduced response to inflammation. LPS was added to study the effects that the introduction of p30 had on HPC-7 capacity to respond to an immune challenge on a transcriptional level. Previous RT-qPCR experiments have shown that after an acute activation (2h) the expression of some acute inflammatory response genes like Ccl5, Cxcl2 or Il12b was not induced as much as in the cells with the long isoform of C/EBPa p42.  

## Samples  

EV 

EV + LPS 

P30 

P30 + LPS 

P42 

P42 + LPS 

 

## Analyses  

### FASTQC 
fastqc.sh

### Align to the mouse genome (mm39) 
hisat2_build_mm39.sh, hisat2_mm39.sh 

### Read counts on GENCODE transcriptome vM27
htseq_count.sh, generate_count_matrix.R  
conversion of ENSEMBL ids to gene names: convert_ENSgene_id_to_gene_name.R  
normalization (FPKM and TPM): calculate_gene_length.R, generate_FPKM_and_TPM.R

### Differential gene expression (DESeq2) and GSEA
Scripts to match ENSEMBL gene IDs and NCBI IDs in mouse gene sets RData files (source: https://bioinf.wehi.edu.au/MSigDB/) and to generate gmt files: create_geneID_conversion_table_for_GSEA.R, generate_gmt_files_from_RData_mouse_gene_sets.R

#### DGEA:  
DESeq2_CEBPa.R


#### GSEA:  
##### Hallmarks:  
generate_bash_run_gsea_Hallmarks.R, run_GseaPreranked_CEBPa_DESeq2.sh  
##### C2, C3, C7:  
generate_bash_run_gsea_c2_c3_c7.R, run_GseaPreranked_CEBPa_DESeq2_c2_c3_c7_gene_sets.sh

#### PCA and Hierarchical clustering 
PCA_and_heatmaps_sample_clustering.R

 

 

 

Other things:  

 

PCA of samples separated between LPS and non-LPS 

Compare induction in different genotypes (LPS vs non LPS) 

Venn diagram of induced genes by p30 and p42 

Bar chart of expression levels before/after induction 

Compare between genotypes 

Plot fold-induction: heatmap, bar charts 

Identify DE pathways before and after inflammation: GSEA, GO 

Identify motif enrichment in DE gene promoters 

Place DE genes into the Bhatt classification 

Identify DE lncRNA (see S Carpenter’s talk on how important are they for inflammation) 

Is there alternative splicing?  

Is there alternative TSS usage? 
