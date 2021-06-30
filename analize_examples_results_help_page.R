dds <- makeExampleDESeqDataSet(n=100,m=18)
dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)
vd <- VisualizeDesign(sampleData = data.frame(condition=dds$condition,genotype=dds$genotype), 
                      designFormula =  ~ genotype + condition + genotype:condition, 
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
app <- ExploreModelMatrix(sampleData = coldata, 
                          designFormula = ~0+ vector + LPS + vector:LPS)
if (interactive()) shiny::runApp(app)
# the condition effect for genotype I (the main effect)

head(results(dds, contrast=c("condition","B","A")))
head(results(dds,contrast = c(0,0,0,1,0,0)))

# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
#BIIIvsAIII:
c(1,0,1,1,0,1) - c(1,0,1,0,0,0) 
BIIIvs_AIII=results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
#note that in our design this is analogous to p42LPSvsp42 

# the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
results(dds, name="genotypeIII.conditionB")

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))
#BIIIvsBII? 
results(dds, contrast = c(0,0,0,0,-1,1))
# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.

# ~~~ Using a grouping variable ~~~

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.

dds$group <- factor(paste0(dds$genotype, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

# the condition effect for genotypeIII
results(dds, contrast=c("group", "IIIB", "IIIA"))
vd <- VisualizeDesign(sampleData = data.frame(condition=dds$condition,genotype=dds$genotype,
                                              group=dds$group), 
                      designFormula =  ~ group, 
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
