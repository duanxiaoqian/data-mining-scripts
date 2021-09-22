# correlation
library(openxlsx)
library(dplyr)
library(DMwR)
library(mice)

setwd("Z:/项目数据/LCMS/2019/LM2019-0944/过程文件")
# 提取数据
PIO<- read.xlsx("送检的卵泡液临床资料！！！.xlsx", rowNames = T)
PO<- read.xlsx("送检的卵泡液临床资料！！！.xlsx", sheet=2, rowNames = T)
P<- read.xlsx("送检的卵泡液临床资料！！！.xlsx", sheet=3,rowNames = T)
c<- read.xlsx("送检的卵泡液临床资料！！！.xlsx", sheet=4,rowNames = T)
# 用knn填充缺失值
mice::md.pattern(PIO)
PIO_knn <- knnImputation(apply(PIO, 2, as.numeric), k=10)
PO_knn <- knnImputation(apply(PO, 2, as.numeric), k=10)
P_knn <- knnImputation(apply(P, 2, as.numeric), k=10)
c_knn <- knnImputation(apply(c, 2, as.numeric), k=10)
merge_C_P_PIO_PO <- rbind(c_knn,P_knn,PIO_knn,PO_knn)

P_C<- read.xlsx("差异代谢物.xlsx")
PIO_PO<- read.xlsx("差异代谢物.xlsx", sheet=2)
PO_P<- read.xlsx("差异代谢物.xlsx", sheet=3)
alldata<- read.xlsx("数据矩阵.xlsx")

# 第一组
one <- alldata[alldata$Metabolites%in%PO_P$Metabolites,][,c(83:102, 43:62)] 
rownames(one) <- PO_P$Metabolites
var <- rbind(PO_knn, P_knn) %>% cbind(., t(one)) %>% cor() %>%
write.xlsx(., "差异代谢物与临床指标的相关分析--PO_P.xlsx", sheetName="差异代谢物与临床指标的相关分析--PO_P", rowNames=T)



# 并集
allmeta <- unique(c(P_C$Metabolites, PIO_PO$Metabolites, PO_P$Metabolites))
all <- alldata[alldata$Metabolites %in% allmeta,]
all <- select(all, -(1:22))
rownames(all) <- allmeta
meta <- t(all)

# merge 各个变量
whole_var <- cbind(merge_C_P_PIO_PO, meta)
whole_cor <- cor(whole_var)

# 保存文件
write.xlsx(whole_cor, "差异代谢物与临床指标的相关分析.xlsx", sheetName="差异代谢物与临床指标的相关分析", rowNames=T)

# all_C <- cbind(select(all, Metabolites),select(all, starts_with("C")))
# all_P <- cbind(select(all, Metabolites),select(all, starts_with("P"))) %>% select(1:21)
# all_PO <- cbind(select(all, Metabolites),select(all, starts_with("PO")))
# all_PIO <- cbind(select(all, Metabolites),select(all, starts_with("PIO")))

# 交集
inmeta <- intersect(intersect(P_C$Metabolites, PIO_PO$Metabolites),PO_P$Metabolites)
inm <- alldata[alldata$Metabolites %in% inmeta,]
inm <- select(inm, -(1:22))
rownames(inm) <- inmeta
meta <- t(inm)
whole_var <- cbind(merge_C_P_PIO_PO, meta)
whole_cor <- cor(whole_var)

save.image()
# 售后1-三组差异代谢物与临床指标分开做相关性分析

#######################################################
# 售后2

setwd("Z:/项目数据/LCMS/2019/LM2019-0944/售后/售后2")

PO<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=1, rowNames = T)
P<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=2,rowNames = T)
c<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=3,rowNames = T)
PO_knn <- knnImputation(apply(PO, 2, as.numeric), k=10)
PO_knn_1 <- PO_knn[,c(22,29:34)]
PO_knn_2 <- PO_knn[, c(1:2,6:7,10)]
P_knn <- knnImputation(apply(P, 2, as.numeric), k=10)
P_knn_1 <- P_knn[,c(22,29:34)]
P_knn_2 <- P_knn[, c(1:2,6:7,10)]
c_knn <- knnImputation(apply(c, 2, as.numeric), k=10)
c_knn_1 <- c_knn[,c(22,29:34)]
c_knn_2 <- c_knn[, c(1:2,6:7,10)]

## diff metabolites
PO_P_MET <- read.xlsx("差异代谢物.xlsx", sheet=2, rowNames = T)
P_C_MET <- read.xlsx("差异代谢物.xlsx", sheet=1, rowNames = T)
all <- read.xlsx("数据矩阵.xlsx", sheet=1, rowNames = F)
PO_P_29 <- all[all$Metabolites%in%PO_P_MET$Metabolites,][,c(83:102,43:62)]
rownames(PO_P_29) <- PO_P_MET$Metabolites
P_C_23 <- all[all$Metabolites%in%P_C_MET$Metabolites,][,c(43:62,23:42)]
rownames(P_C_23) <- P_C_MET$Metabolites
## merge data and cor P
library(psych)
#1
PO_P_1 <- cbind(t(PO_P_29), rbind(PO_knn_1, P_knn_1))
cor1 <- corr.test(PO_P_1, adjust = "none")

RMETA2::savexlsx1(cor1$r, "相关性和pvalue.xlsx", sheet="PO_P_1_cor")
RMETA2::savexlsx1(cor1$p, "相关性和pvalue.xlsx", sheet="PO_P_1_pvalue")

#2
P_C_1 <- cbind(t(P_C_23), rbind(P_knn_1, c_knn_1))
cor2 <- corr.test(P_C_1, adjust = "none")

RMETA2::savexlsx1(cor2$r, "相关性和pvalue.xlsx", sheet="P_C_1_cor")
RMETA2::savexlsx1(cor2$p, "相关性和pvalue.xlsx", sheet="P_C_1_pvalue")

#3
PO_P_2 <- cbind(t(PO_P_29), rbind(PO_knn_2, P_knn_2))
cor3 <- corr.test(PO_P_2, adjust = "none")

RMETA2::savexlsx1(cor3$r, "相关性和pvalue.xlsx", sheet="PO_P_2_cor")
RMETA2::savexlsx1(cor3$p, "相关性和pvalue.xlsx", sheet="PO_P_2_pvalue")

#4
P_C_2 <- cbind(t(P_C_23), rbind(P_knn_2, c_knn_2))
cor4 <- corr.test(P_C_2, adjust = "none")

RMETA2::savexlsx1(cor4$r, "相关性和pvalue.xlsx", sheet="P_C_2_cor")
RMETA2::savexlsx1(cor4$p, "相关性和pvalue.xlsx", sheet="P_C_2_pvalue")

##PCA（all）(剔除PIO组)
setwd("Z:/项目数据/LCMS/2019/LM2019-0944/售后/售后2/要求1")

RMETA2::createregistration()                 #创建项目登记单
data<-RMETA2::readregistration()     #提取项目登记单信息
# data<-RMETA2::addregistration(data,name="项目登记单1.xlsx")
data<-RMETA2::autodataprocess(data)  #获取数据并处理,rsd=NA,rsd=NA,missvalue = NA
#1.数据获取    data<-RMETA2::datafetch(data)
#              data$data$rawdata$data<-RMETA2::readxlsx(name="数据矩阵.xlsx",sheet=1)  #直接有数据矩阵
#另外获取模式  data<-RMETA2::dataregistration(name="数据矩阵.xlsx",sheet="数据矩阵",type="LCMS",dif = ,species = ) 
#2.预处理      data<-RMETA2::preprocess(data)
#3.缺省值筛选  data<-RMETA2::automissvalue(data)
#4.rsd筛选     data<-RMETA2::autorsdfilter(data)
#5.0值处理     data<-RMETA2::zeroprocess(data)
# RMETA2::getdata(data) #矩阵获取 
#最终矩阵获取  RMETA2::getfinaldata(data)
#选择矩阵获取  RMETA2::getselectdata(data,which="data")
# data[["info"]][["mulstatistics"]][["shape"]][[1]]<-c(21,22,23,25,24)
# data[["info"]][["mulstatistics"]][["fill"]][[1]]<-c("green4","blue3", "gold","darkviolet","firebrick")
data<-RMETA2::allstatistics(data)          #pca(all)运算
RMETA2::getallstatisticsmap(data,showname=F)  #获取有名字的pca(all)图
#RMETA2::getstatisticsmap(data,showname=F)
 
## heatmap--质控
RMETA2::createlist(data)

##2 3group heatmap

RMETA2::savexlsx3(P_C_23, "heatmap.xlsx", "P vs C")
RMETA2::savexlsx3(PO_P_29, "heatmap.xlsx", "PO vs p")

union <- unique(c(rownames(P_C_23), rownames(PO_P_29)))
union1 <- all[all$Metabolites%in%union,][, c(5, 83:102,23:62)]
RMETA2::savexlsx1(union1, "heatmap.xlsx", "P_PO_C_union")
intersect <- multi_intersect(list(rownames(P_C_23), rownames(PO_P_29)))
P_PO_C_intersect <- all[all$Metabolites%in%intersect,][, c(5, 83:102,23:62)]
RMETA2::savexlsx1(P_PO_C_intersect, "heatmap.xlsx", "P_PO_C_intersect")
RMETA2::auto_alldraw(name = 'heatmap.xlsx',
                     drawwhat = "heatmap",row = F,col =F) 

#60个样的PCA
PCA_union <- all[all$Metabolites%in%union,]
PCA_intersect <- all[all$Metabolites%in%intersect,]
RMETA2::savexlsx1(PCA_union, "PCA.xlsx", "PCA_union")
RMETA2::savexlsx1(PCA_intersect, "PCA.xlsx", "PCA_intersect")

data<-RMETA2::readregistration()     #提取项目登记单信息
# data<-RMETA2::addregistration(data,name="项目登记单1.xlsx")
data<-RMETA2::autodataprocess(data)  #获取数据并处理,rsd=NA,rsd=NA,missvalue = NA
#1.数据获取    data<-RMETA2::datafetch(data)
#              data$data$rawdata$data<-RMETA2::readxlsx(name="PCA.xlsx",sheet=2)  #直接有数据矩阵
#另外获取模式  data<-RMETA2::dataregistration(name="数据矩阵.xlsx",sheet="数据矩阵",type="LCMS",dif = ,species = ) 
#2.预处理      data<-RMETA2::preprocess(data)
#3.缺省值筛选  data<-RMETA2::automissvalue(data)
#4.rsd筛选     data<-RMETA2::autorsdfilter(data)
#5.0值处理     data<-RMETA2::zeroprocess(data)
# RMETA2::getdata(data) #矩阵获取 
#最终矩阵获取  RMETA2::getfinaldata(data)
#选择矩阵获取  RMETA2::getselectdata(data,which="data")
# data[["info"]][["mulstatistics"]][["shape"]][[1]]<-c(21,22,23,25,24)
# data[["info"]][["mulstatistics"]][["fill"]][[1]]<-c("green4","blue3", "gold","darkviolet","firebrick")
data<-RMETA2::allstatistics(data)          #pca(all)运算
RMETA2::getallstatisticsmap(data,showname=F)


#2.kegg
setwd("Z:/项目数据/LCMS/2019/LM2019-0944/售后/售后2/要求2")
PO_P_KEGG <- all[all$Metabolites%in%rownames(PO_P_29),]
P_C_KEGG <- all[all$Metabolites%in%rownames(P_C_23),]
RMETA2::savexlsx1(PO_P_KEGG, "KEGG_enrichment.xlsx", "PO_P_KEGG")
RMETA2::savexlsx1(P_C_KEGG, "KEGG_enrichment.xlsx", "P_C_KEGG")

RMETA2::auto_allweb(name="KEGG_enrichment.xlsx",KEGG="local",rich="local",
                    from="LCMS",needlist=c("Metabolites","Compound.ID"),needgroup=NULL,species="hsa")


# 售后3
## heatmap
setwd("Z:/项目数据/LCMS/2019/LM2019-0944/售后/售后3")
P_C <- RMETA2::readxlsx("差异代谢物1.xlsx", sheet=1) 
P_C_40 <- P_C[,c(5, 22:61)]
P_C_avg <- P_C[,c(5, 20:21)]
colnames(P_C_avg) <- c("Metabolites", "P", "C")
RMETA2::savexlsx1(P_C_40, "heatmap.xlsx", "P_C_40")
RMETA2::savexlsx1(P_C_avg, "heatmap.xlsx", "P_C_avg")

PO_P <- RMETA2::readxlsx("差异代谢物1.xlsx", sheet=2) 
PO_P_40 <- PO_P[, c(5, 22:61)]
PO_P_avg <- PO_P[, c(5, 20:21)]
colnames(PO_P_avg) <- c("Metabolites", "PO", "P")
RMETA2::savexlsx1(PO_P_40, "heatmap.xlsx", "PO_P_40")
RMETA2::savexlsx1(PO_P_avg, "heatmap.xlsx", "PO_P_avg")

RMETA2::auto_alldraw(name = 'heatmap.xlsx',
                     drawwhat = "heatmap",row = F,col =F)
# correlation 
## 临床数据处理
library(openxlsx)
library(psych)
PO<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=1, rowNames = T)
P<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=2,rowNames = T)
c<- read.xlsx("附件2：临床数据--LM2019-0944.xlsx", sheet=3,rowNames = T)
PO_knn <- knnImputation(apply(PO, 2, as.numeric), k=10)[,c(1:2, 6:7, 10, 22, 29:34)]
P_knn <- knnImputation(apply(P, 2, as.numeric), k=10)[,c(1:2, 6:7, 10, 22, 29:34)]
c_knn <- knnImputation(apply(c, 2, as.numeric), k=10)[,c(1:2, 6:7, 10, 22, 29:34)]

## P_1
P_MET <- P_C[, c(22:41)]
rownames(P_MET) <- P_C$Metabolites
P_COR_1 <- cbind(t(P_MET),P_knn) %>% corr.test(.,adjust = "none")
RMETA2::savexlsx3(P_COR_1$r, "相关性和pvalue.xlsx", sheet="P_COR_1")
RMETA2::savexlsx3(P_COR_1$p, "相关性和pvalue.xlsx", sheet="P_pvalue_1")

## C_1
C_MET <- P_C[, c(42:61)]
rownames(C_MET) <- P_C$Metabolites
C_COR_1 <- cbind(t(C_MET),c_knn) %>% corr.test(.,adjust = "none")
RMETA2::savexlsx3(C_COR_1$r, "相关性和pvalue.xlsx", sheet="C_COR_1")
RMETA2::savexlsx3(C_COR_1$p, "相关性和pvalue.xlsx", sheet="C_pvalue_1")

## PO_1
PO_MET <- PO_P[, c(22:41)]
rownames(PO_MET) <- PO_P$Metabolites
PO_COR_1 <- cbind(t(PO_MET),PO_knn) %>% corr.test(.,adjust = "none")
RMETA2::savexlsx3(PO_COR_1$r, "相关性和pvalue.xlsx", sheet="PO_COR_1")
RMETA2::savexlsx3(PO_COR_1$p, "相关性和pvalue.xlsx", sheet="PO_pvalue_1")

## P_2
P_MET_2 <- PO_P[, c(42:61)]
rownames(P_MET_2) <- PO_P$Metabolites
P_COR_2 <- cbind(t(P_MET_2),P_knn) %>% corr.test(.,adjust = "none")
RMETA2::savexlsx3(P_COR_2$r, "相关性和pvalue.xlsx", sheet="P_COR_2")
RMETA2::savexlsx3(P_COR_2$p, "相关性和pvalue.xlsx", sheet="P_pvalue_2")
