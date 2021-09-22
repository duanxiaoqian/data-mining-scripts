# after-sevice for 1475 -1 -2 -3
setwd("Z:/项目数据/LCMS/2019/LM2019-1475-售后/售后1")
library(dplyr)
# heatmap, corplot
ALL <- RMETA2::readxlsx("数据矩阵-肝脏.xlsx",sheet = 1)
MOD_CON <- RMETA2::readxlsx("筛选差异代谢物-肝脏 (1).xlsx",sheet = 1)
MOD_QGD <- RMETA2::readxlsx("筛选差异代谢物-肝脏 (1).xlsx",sheet = 2)
MOD_HGD <- RMETA2::readxlsx("筛选差异代谢物-肝脏 (1).xlsx",sheet = 3)
MOD_CP <- RMETA2::readxlsx("筛选差异代谢物-肝脏 (1).xlsx",sheet = 4)

MOD_CON_1 <- cbind(dplyr::select(MOD_CON, Metabolites),
                   select(MOD_CON, starts_with("MOD")), 
                   select(MOD_CON, starts_with("CON")))
RMETA2::savexlsx1(MOD_CON_1, name = "heatmap.xlsx", "MOD_CON")
MOD_QGD_1 <- cbind(select(MOD_QGD, Metabolites),
                   select(MOD_QGD, starts_with("MOD")), 
                   select(MOD_QGD, starts_with("QGD")))
RMETA2::savexlsx1(MOD_QGD_1, name = "heatmap.xlsx", "MOD_QGD")
MOD_HGD_1 <- cbind(select(MOD_HGD, Metabolites),
                   select(MOD_HGD, starts_with("MOD")), 
                   select(MOD_HGD, starts_with("HGD")))
RMETA2::savexlsx1(MOD_HGD_1, name = "heatmap.xlsx", "MOD_HGD")
MOD_CP_1 <- cbind(select(MOD_CP, Metabolites),
                  select(MOD_CP, starts_with("MOD")), 
                  select(MOD_CP, starts_with("CP")))
RMETA2::savexlsx1(MOD_CP_1, name = "heatmap.xlsx", "MOD_CP")

all_me <- unique(c(MOD_CON_1$Metabolites, MOD_QGD_1$Metabolites,
            MOD_HGD_1$Metabolites, MOD_CP_1$Metabolites))
ALL$QGD10
all_1 <- cbind(select(ALL, Metabolites),
               select(ALL, CON6:QGD10))
all_2 <- all_1[,unique(c(colnames(MOD_CON_1), colnames(MOD_QGD_1), colnames(MOD_HGD_1), colnames(MOD_CP_1)))]

comm_5group <- all_2[all_2$Metabolites%in%all_me,]
RMETA2::savexlsx1(comm_5group, name = "heatmap.xlsx", "ALL_5Group")

# inter_me <- intersect(intersect(MOD_CON_1$Metabolites, MOD_QGD_1$Metabolites),
#             intersect(MOD_HGD_1$Metabolites, MOD_CP_1$Metabolites))
# inter_all <- cbind(MOD_CON_1[MOD_CON_1$Metabolites%in%inter_me, ],
#                    MOD_QGD_1[MOD_QGD_1$Metabolites%in%inter_me, ] %>% select(., starts_with("QGD")),
#                    MOD_HGD_1[MOD_HGD_1$Metabolites%in%inter_me, ] %>% select(., starts_with("HGD")),
#                    MOD_CP_1[MOD_CP_1$Metabolites%in%inter_me, ] %>% select(., starts_with("CP")))
# RMETA2::savexlsx1(inter_all, name = "heatmap.xlsx", "ALL_5Group")

RMETA2::auto_alldraw(name="heatmap.xlsx",drawwhat="heatmap",row=F,col=F,datasave = F) #给热图矩阵画热图
RMETA2::auto_alldraw(name="heatmap-血清.xlsx", drawwhat="correlation",datasave = F)

# kegg analysis

RMETA2::auto_allweb(name="筛选差异代谢物-血清 (1).xlsx",KEGG="local",rich="local",from="LCMS",
                    needlist=c("Metabolites","Compound ID","FC"),needgroup=NULL,species="rno")

# loading  plot, S plot, VIP plot
RMETA2::getloadingsmap(data,showname = F,mode = "OPLS", wbsave = T) # data save and loading plot
RMETA2::getsplotmap(data,mode = "OPLS",showname=T,wbsave=T,showpar=T,type=c("jpg","pdf"))
# get vip data
library(RMETA2)
RMETA2::savevip(data) 
library(ggplot2)
library(openxlsx)

sheetname <- openxlsx::getSheetNames("VIP.xlsx")
for (i in 1:length(sheetname)){
  data1 <- RMETA2::readxlsx("VIP.xlsx", sheet = i) %>% arrange(., desc(VIP))
  data1 <- data1[1:300,]
  p <- ggplot(data=data1, aes(x=1:length(data1$VIP), y=VIP)) + 
  geom_col(fill="#339900")+
  xlab("Num") + ylab("VIP")+ggtitle(paste0(sheetname[i], " VIP Plot"))+
  theme_bw()+
  theme(panel.grid=element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  scale_x_continuous(breaks = seq(0,300,20))+
  scale_y_continuous(breaks = seq(0,50,5))
  pdf(paste0(sheetname[i], "-VIP_Plot.pdf"),width=15, height=8)
  print(p)
  dev.off()
  jpeg(paste0(sheetname[i], "-VIP_Plot.jpg"),width=900, height=500)
  print(p)
  dev.off()
  print("VIP Plot drawing done")
}
  

  
dev.off()


