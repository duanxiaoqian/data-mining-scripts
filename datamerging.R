
setwd("Z:/项目数据/LCMS/2019/LM2019-0628/售后1")
library(openxlsx)
library(dplyr)

data_D_DH <- openxlsx::read.xlsx("差异代谢物.xlsx", sheet="CaseDH_ConD")
data_S_SH <- openxlsx::read.xlsx("差异代谢物.xlsx", sheet="CaseSH_ConS")
data_all <- openxlsx::read.xlsx("数据矩阵.xlsx", sheet=1)

data_common <- merge(data_D_DH, data_S_SH, by = "Metabolites",incomparables = NA)
data <- data_all[data_all[,"Metabolites"] %in% data_common[,"Metabolites"], ]
# data_merged <- cbind(data %>% {c(select(.,Metabolites),
#                                 select(.,starts_with("D")),
#                                 select(.,starts_with("S-")),
#                                 select(.,starts_with("SH")))}) # reture list type
data_merged <- cbind(data %>% select(Metabolites), 
                     data %>% select(starts_with("D")),   
                     data %>% select(starts_with("S-")),
                     data %>% select(starts_with("SH")))
RMETA2::savexlsx1(data, "共有差异代谢物数据矩阵.xlsx", sheet="共有差异代谢物数据矩阵")
RMETA2::savexlsx1(data_merged, "共有差异代谢物.xlsx", sheet="共有差异代谢物")
# 画图
RMETA2::auto_alldraw(name="共有差异代谢物.xlsx",needgroup=NULL,drawwhat="heatmap",row=T,col=F,datasave = F,range = NULL) #给热图矩阵画热图

# 共有代谢物，同上，同下
data_D_DH_UP <- data_D_DH[data_D_DH$`log2(FC)`>0,]
data_D_DH_down <- data_D_DH[data_D_DH$`log2(FC)`<0,]
data_S_SH_up <- data_S_SH[data_S_SH$`log2(FC)`>0,]
data_S_SH_down <- data_S_SH[data_S_SH$`log2(FC)`<0,]
up <- data_D_DH_UP[data_D_DH_UP$Metabolites %in% data_S_SH_up$Metabolites,]
down <- data_D_DH_down[data_D_DH_down$Metabolites %in% data_S_SH_down$Metabolites,]

upanddown <- data_merged[data_merged$Metabolites %in% up$Metabolites | data_merged$Metabolites %in% down$Metabolites,]
RMETA2::savexlsx1(up, "共有上调差异代谢物.xlsx", sheet="共有上调差异代谢物")
RMETA2::savexlsx1(down, "共有下调差异代谢物.xlsx", sheet="共有下调差异代谢物")
RMETA2::savexlsx1(upanddown, "共有上调+下调差异代谢物.xlsx", sheet="共有上调+下调差异代谢物")
# kegg
RMETA2::auto_allweb(name="共有上调差异代谢物.xlsx",KEGG="local",rich="local",from="LCMS",needlist=c("Metabolites","Compound.ID","log2(FC)"),needgroup=NULL,species="mtr")
RMETA2::auto_allweb(name="共有下调差异代谢物.xlsx",KEGG="local",rich="local",from="LCMS",needlist=c("Metabolites","Compound.ID","log2(FC)"),needgroup=NULL,species="mtr")


# corrplot
RMETA2::auto_alldraw(name="共有差异代谢物.xlsx", drawwhat="correlation", datasave = T)
RMETA2::auto_alldraw(name="共有上调+下调差异代谢物.xlsx", drawwhat="correlation", datasave = T)
