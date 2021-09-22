setwd("Z:/项目数据/LCMS/2019/LM2019-0720/售后/售后1")
library(openxlsx)
library(dplyr)

data <- openxlsx::read.xlsx("venn-L_N&H_N&L_H.xlsx", sheet="并集") 

data_heatmap <- cbind(data %>% select(Metabolites), data %>% select(starts_with("average")))
colnames(data_heatmap) <- c("Metabolites", "N", "L","H")
RMETA2::savexlsx1(data_heatmap, "并集-heatmap.xlsx", sheet="并集-heatmap")

# 画图
RMETA2::auto_alldraw(name="并集-heatmap.xlsx",needgroup=NULL,drawwhat="heatmap",row=T,col=F,datasave = F,range = NULL)

# 时序分析
RMETA2::dSTEM(filename = "并集-heatmap.xlsx",sheet = 1, cluster = 15)
