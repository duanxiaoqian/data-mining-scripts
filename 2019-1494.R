# 售后LM2019-1494
setwd("Z:/项目数据/LCMS/2020/LM2019-1494/售后/售后1")

library(dplyr)
#veen
d3 <- RMETA2::readxlsx("veen.xlsx", sheet = 1) 
d3_filtered <- d3[apply(d3[,-1],1, FUN = function(x){return(sum(x==0)<3)}),] # 比较时，必须全为numeric，不然全是false

d6<- RMETA2::readxlsx("veen.xlsx", sheet = 2)
d6_filtered <- d6[apply(d6[,-1], 1, FUN = function(x){return(sum(x==0)<3)}),]

RMETA2::savexlsx1(d3_filtered, "data_filtered.xlsx",sheet="d3_filtered")
RMETA2::savexlsx1(d6_filtered, "data_filtered.xlsx",sheet="d6_filtered")

RMETA2::auto_alldraw(name="data_filtered.xlsx",needgroup=list(c(1,2)),drawwhat="venny",needlist="Metabolites")
RMETA2::auto_alldraw(name = 'heatmap.xlsx',
                     drawwhat = "heatmap",row = F,col =F) # 画行 列都聚类的图

