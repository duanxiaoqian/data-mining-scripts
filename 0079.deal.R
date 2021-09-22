# 0079售后3
library(RMETA2)
library(dplyr)

setwd("Z:/项目数据/脂质组学/2019/LM2019-0079-小鼠血清/售后/售后4")
alldata_VIP <- readxlsx("非差异代谢物-VIP&P.xlsx")
alldata_FC <- readxlsx("非差异代谢物-FC&P.xlsx")

FC <- readxlsx("差异代谢物-FC&P.xlsx")
VIP <- readxlsx("差异代谢物-VIP&P.xlsx")

noFC <- alldata_FC[!(alldata_FC$Metabolites%in%FC$Metabolites),] %>% 
        select(., -c(1:5,7,9:18)) %>% dplyr::arrange(., Class, .by_group=T)
noVIP<- alldata_VIP[!(alldata_VIP$Metabolites%in%VIP$Metabolites),] %>% 
  select(., -c(1:5,7,9:19)) %>% dplyr::arrange(., Class, .by_group=T)

length(unique(noFC$Class))
length(unique(noVIP$Class))

RMETA2::savexlsx1(noFC,"所有非差异代谢物-FC&P.xlsx","所有非差异代谢物-FC&P")
RMETA2::savexlsx1(noVIP,"所有非差异代谢物-VIP&P.xlsx","所有非差异代谢物-VIP&P")

ele_number <- function(Vector=Vector, elements=elements,...){
  number <- c()
  for (i in 1:length(elements)){
    number[i] <- sum(Vector%in%elements[i])
  }
  pd <- data.frame(elements=elements,
                   number=number)
  return(pd)
}

noFC_sort <- ele_number(noFC$Class, unique(noFC$Class)) %>% arrange(.,desc(number))
j<-0
for (i in noFC_sort$elements){
  noFC_ele <- rbind(filter(noFC, noFC$Class==i)[,-2], c("Colgroup", rep("WT",6), rep("CKO",6)))
  RMETA2::savexlsx1(noFC_ele, name = "非差异代谢物-FC&P-Class分类.xlsx",paste0("非差异代谢物-FC&P-",i))
  j=j+1
  print(j)
}
  RMETA2::auto_alldraw(name="非差异代谢物-FC&P-Class分类.xlsx",
                       drawwhat="heatmap",
                       show_rownames = T,
                       # rowgroup = 1,
                       colgroup =1)

nchar()
