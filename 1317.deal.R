# 1317 saleback

library(RMETA2)


setwd("Z:/home/段光前/20191210-LM2019-1317/售后/售后1")
# volcane plot
RMETA2::auto_alldraw(name="volcano-metabolites.xlsx",drawwhat="volcano",type=c("html","pdf","jpg"))
data <- RMETA2::readxlsx("volcano-metabolites.xlsx",sheet = 1, rowNames = T)
auto_volcano(data = data,
             name = "SA_Con",
             type = c("jpg","pdf","html"))

# heatmap
RMETA2::auto_alldraw(name = './heatmap.xlsx',
                     drawwhat = "heatmap",row = F,col =F ,type=c("jpg","pdf"), 
                     color = colorRampPalette(c("#0066CC","#FFFFFF", "#FF3300"))(10000))

dev.off()

# correlation  plot
data_me1 <- RMETA2::readxlsx("差异代谢物.xlsx", sheet=1)
data_me2 <- RMETA2::readxlsx("差异代谢物.xlsx", sheet=2)
data_cellfa <- RMETA2::readxlsx("细胞因子.xlsx", sheet=1,rowNames = T)

data_me2 <- select(data_me2, one_of("Metabolites", "SA1",
                                                   "SA3","SA6","SA11","SA17",
                                                   "SA21","SA24","SA27","SA32","SA34")) %>% t()
data_cellfa <-  data_cellfa[c("SA1",
                              "SA3","SA6","SA11","SA17",
                              "SA21","SA24","SA27","SA32","SA34"),]                
colnames(data_me2) <- data_me2[1,]
data_me2 <- data_me2[-1,]
data_me2_merge <- cbind(data_cellfa,data_me2) %>% 
  apply(., 2, as.numeric) %>% psych::corr.test() 
data_re <- data_me2_merge[["r"]] %>% RMETA2::savexlsx1(., name = "相关性矩阵.xlsx", sheet="SA-1_SA-2")
data_me2_merge[["p"]] %>% RMETA2::savexlsx1(., name = "显著性矩阵.xlsx", sheet="SA-1_SA-2")











library(magrittr)
library(dplyr)
library(openxlsx)
data_me1 <- cbind(select(data_me1,Metabolites),
                  select(data_me1,starts_with("Con")),
                  select(data_me1, starts_with("SA"))) %>% t() 
colnames(data_me1) <- data_me1[1,]
data_me1 <- data_me1[-1,]
data_me1_merge <- cbind(data_cellfa,data_me1) %>% 
  apply(., 2, as.numeric) %>% psych::corr.test() 
data_re <- data_me1_merge[["r"]]  RMETA2::savexlsx1(., name = "相关型矩阵.xlsx", sheet="SA_Con")
data_me1_merge[["p"]] %>% RMETA2::savexlsx1(., name = "显著性矩阵.xlsx", sheet="SA_Con")
data <- RMETA2::readxlsx("相关型矩阵.xlsx", rowNames=T)
library(ComplexHeatmap)
# h <- ComplexHeatmap::Heatmap(data,
#                    color=colorRampPalette(c("#0099FF","#FFFFFF", "#FF6699"))(1000),
#                    cluster_columns  = FALSE,
#                    
#                    )
library(pheatmap)
 P <- pheatmap(data, color = colorRampPalette(c("#0099FF","#FFFFFF", "#FF6699"))(1000),
         fontsize_row=2,
         cluster_cols = F)
pdf("correlation.pdf",width = 10,height = 12)
print(P)
dev.off()
RMETA2::auto_alldraw(name = './相关性矩阵.xlsx',
                     drawwhat = "heatmap",row = F,col =F ,type=c("jpg","pdf"), 
                     color = colorRampPalette(c("#0099FF","#FFFFFF", "#FF6699"))(1000),
                     scale = "none",
                     cluster_rows = T)

draw(h)
dev.off()

# veen
RMETA2::auto_alldraw(name="差异代谢物.xlsx",needgroup=list(c(1:2)),drawwhat="venny",needlist="Metabolites")


# VIP
RMETA2::auto_alldraw(name="差异代谢物.xlsx",needgroup=NULL,drawwhat="vip",needlist=c("Metabolites", "VIP", "SA-1","SA-2"))
