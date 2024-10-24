p42LPS_vs_p30LPS=results(dds,contrast = c(0,-1,1,0,-1,1))
head(p42LPS_vs_p30LPS)
summary(p42LPS_vs_p30LPS)

p42_vs_p30 <- results(dds,contrast = c(0,-1,1,0,0,0))
summary(p42_vs_p30)

p42_vs_p30_alt <- results(dds, contrast = c("vector","p42","p30"))
summary(p42_vs_p30_alt)

p30_vs_EV=results(dds,contrast = c(-1,1,0,0,0,0))
summary(p30_vs_EV)
p30_vs_EV_alt <- results(dds,contrast = c("vector","p30","EV"))
summary(p30_vs_EV)

p42_vs_EV=results(dds,contrast = c(-1,0,1,0,0,0))
summary(p42_vs_EV)
p42_vs_EV_alt <- results(dds,contrast = c("vector","p42","EV"))
summary(p42_vs_EV)

#P30_LPS vs EV_LPS 
p30LPS_vs_EVLPS <- results(dds,c(-1,1,0,0,1,0))
summary(p30LPS_vs_EVLPS)
p30LPS_vs_EVLPS_alt <- results(dds,contrast = list())

#p42_LPS vs EV_LPS 
p42LPS_vs_EVLPS <- results(dds,c(-1,0,1,0,0,1))
summary(p42LPS_vs_EVLPS)
p42LPS_vs_EVLPS_alt <- results(dds,contrast = list(c("vectorp42" ,"vectorp42.LPSTRUE"),"vectorEV"))
summary(p42LPS_vs_EVLPS_alt)


#EVLPS vs EV
EVLPS_vs_EV <- results(dds,c(0,0,0,1,0,0))
summary(EVLPS_vs_EV)
EVLPS_vs_EV_alt <- results(dds,contrast = c("LPS","TRUE","FALSE"))
summary(EVLPS_vs_EV_alt)
EVLPS_vs_EV_alt2 <- results(dds,name = c("LPSTRUE"))
summary(EVLPS_vs_EV_alt2)

#P30_LPS vs P30
p30LPS_vs_p30 <- results(dds,contrast = c(0,0,0,1,1,0))
summary(p30LPS_vs_p30)
#another way of specifying this contrast is:

p30LPS_vs_p30_alt <- results(dds,contrast = list(c("vectorp30.LPSTRUE","LPSTRUE")))
summary(p30LPS_vs_p30_alt)
list(c("vectorp30.LPSTRUE","LPSTRUE"))
#note this: if the list is length 1, a second element is added which is the empty character vector, character()

#If I mistakenly forget to include the c() then this is the result
p30LPS_vs_p30_wrong <- results(dds,contrast = list("vectorp30.LPSTRUE","LPSTRUE"))
summary(p30LPS_vs_p30_wrong)
list("vectorp30.LPSTRUE","LPSTRUE")

#so it is telling the function to consider "vectorp30.LPSTRUE" as numerator
#and "LPSTRUE" as denominator

#this should be equivalent to:
test <- results(dds,contrast = c(0,0,0,-1,1,0))
summary(test)
#it is indeed the same

#P42_LPS vs p42 
p42LPS_vs_p42 <- results(dds,contrast=c(0,0,0,1,0,1))
summary(p42LPS_vs_p42)


#taken from ?results:
# ~~~ Using a grouping variable ~~~

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.
#dds$group <- factor(paste0(dds$genotype, dds$condition))
#design(dds) <- ~ group
ddsMat$group=factor(paste0(ddsMat$vector,ddsMat$LPS))
design(ddsMat) <- ~group

dds <- DESeq(ddsMat)
resultsNames(dds)
p42LPS_vs_p42_groupdesign <- results(dds, contrast = c("group","p42TRUE","p42FALSE"))
summary(p42LPS_vs_p42_groupdesign)
summary(p42LPS_vs_p42)
#same result!!
BiocManager::install("ExploreModelMatrix")
library(ExploreModelMatrix)
vd <- VisualizeDesign(sampleData = data.frame(group=dds$group), 
                      designFormula =  ~ group, 
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
#p30vsEV
p30_vs_EV_gd <- results(dds,contrast = c("group", "p30FALSE","EVFALSE"))
summary(p30_vs_EV)
summary(p30_vs_EV_gd)



# P42 vs EV 
# 
p42_vs_EV_gd <- results(dds,contrast = c("group", "p42FALSE","EVFALSE"))
summary(p42_vs_EV_gd)
summary(p42_vs_EV)
# P30 vs p42 
# 
p42_vs_p30_gd <- results(dds,contrast = c("group", "p42FALSE","p30FALSE"))
summary(p42_vs_p30)
summary(p42_vs_p30_gd)

# P30_LPS vs EV_LPS 
# 
p30LPS_vs_EVLPS_gd <- results(dds,contrast = c("group", "p30TRUE","EVTRUE"))
summary(p30LPS_vs_EVLPS)
summary(p30LPS_vs_EVLPS_gd)
# P42_LPS vs EV_LPS 
# 
p42LPS_vs_EVLPS_gd <- results(dds,contrast = c("group", "p42TRUE","EVTRUE"))
summary(p42LPS_vs_EVLPS)
summary(p42LPS_vs_EVLPS_gd)
# P30_LPS vs P42_LPS 
# 
p42LPS_vs_p30LPS_gd <- results(dds,contrast = c("group", "p42TRUE","p30TRUE"))
summary(p42LPS_vs_p30LPS)
summary(p42LPS_vs_p30LPS_gd)
# EV_LPS vs EV 
# 
EVLPS_vs_EV_gd <- results(dds,contrast = c("group", "EVTRUE","EVFALSE"))
summary(EVLPS_vs_EV_gd)
summary(EVLPS_vs_EV)
# P30_LPS vs P30 
# 
p30LPS_vs_p30_gd <- results(dds,contrast = c("group", "p30TRUE","p30FALSE"))
summary(p30LPS_vs_p30)
summary(p30LPS_vs_p30_gd)
# P42_LPS vs p42 
# 
p42LPS_vs_p42_gd <- results(dds,contrast = c("group","p42TRUE","p42FALSE"))
summary(p42LPS_vs_p42)
summary(p42LPS_vs_p42_gd)


#so, for the type of comparisons we want to perform, it is the same to
#indicate the design as groups. I think the potential problem here is that
#we do not have the interaction term in this case

#remember:
# ~~~ Using a grouping variable ~~~
#"# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.
#"

#and that is exactly what we want to compare