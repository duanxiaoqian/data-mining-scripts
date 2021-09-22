# LM2019-1650-售后1
setwd("Z:/")

# read data

BFB_BFW <- RMETA2::readxlsx("差异代谢物.xlsx",sheet=1)
BFB_NBF <- RMETA2::readxlsx("差异代谢物.xlsx",sheet=2)
BFW_NBF <- RMETA2::readxlsx("差异代谢物.xlsx",sheet=3)
all <- RMETA2::readxlsx("数据矩阵.xlsx",sheet=1)

# union
union <- unique(c(BFB_BFW$Metabolites, BFB_NBF$Metabolites, BFW_NBF$Metabolites))

# intersect
intersect <- multi_intersect(list(BFB_BFW$Metabolites, BFB_NBF$Metabolites, BFW_NBF$Metabolites))

# heatmap
heatmap_union <- all[all$Metabolites%in%union, ][, c(5, 18:26)]
heatmap_intersect <- all[all$Metabolites%in%intersect, ][, c(5, 18:26)]
RMETA2::savexlsx1(heatmap_union, "heatmap.xlsx", sheet = "heatmap_union")
RMETA2::savexlsx1(heatmap_intersect, "heatmap.xlsx", sheet = "heatmap_intersect")

RMETA2::auto_alldraw(name = 'heatmap.xlsx',
                     drawwhat = "heatmap",row = F,col =F) 
# CORRELATION
RMETA2::auto_alldraw(name = 'heatmap.xlsx',
                     drawwhat = "correlation", datasave = T) 

# kegg
kegg_union <- all[all$Metabolites%in%union, ]
kegg_intersect <- all[all$Metabolites%in%intersect, ]
RMETA2::savexlsx1(kegg_union, "kegg.xlsx", sheet = "kegg_union")
RMETA2::savexlsx1(kegg_intersect, "kegg.xlsx", sheet = "kegg_intersect")
RMETA2::auto_allweb(name="kegg.xlsx",KEGG="local",rich="local",from="LCMS",
  needlist=c("Metabolites","Compound ID"),needgroup=NULL,species="mdm")

