setwd("Z:/项目数据/LCMS/2019/LM2019-0637/售后2/1.1 单变量统计分析-paired T test")
####1
# heatmap data extract
library(openxlsx)
library(dplyr)

# heatmap
data_8 <- openxlsx::read.xlsx("经likelihood ratio test调整的混合效应识别差异代谢物数据矩阵—17种.xlsx")
data_8heatmap <- cbind(select(data_8,Metabolites), select(data_8,starts_with("A1")), select(data_8,starts_with("A2")))
RMETA2::savexlsx1(data_8heatmap,name = "heatmap.xlsx",sheet = "heatmap")
RMETA2::auto_alldraw(name="heatmap.xlsx",needgroup=NULL,drawwhat="heatmap",row=T,col=F,datasave = F,range = NULL)

# corrplot
RMETA2::auto_alldraw(name="heatmap.xlsx",drawwhat="correlation",datasave = T)

# kegg
# 添加log2FC
data_8[,"log2(FC)"] <- log2(rowMeans(select(data_8,starts_with("A2")))/rowMeans(select(data_8,starts_with("A1"))))
colnames(data_8)
RMETA2::savexlsx1(data_8,name = "经likelihood ratio test调整的混合效应识别差异代谢物数据矩阵—17种+FC.xlsx",sheet = "差异代谢物数据矩阵")
RMETA2::auto_allweb(name="经likelihood ratio test调整的混合效应识别差异代谢物数据矩阵—17种+FC.xlsx",KEGG="local",rich="local",from="LCMS",
                    needlist=c("Metabolites","Compound.ID","log2(FC)"),needgroup=NULL,species="hsa")


#######2
# P<0.05差异代谢物
data_0.05 <- openxlsx::read.xlsx("P0.05差异代谢物.xlsx")
data_301 <- openxlsx::read.xlsx("线性混合效应模型数据矩阵0923.xlsx")
data_301heatmap <- cbind(select(data_301[data_301$Metabolites %in% data_0.05$Metabolites,],Metabolites), 
                         select(data_301[data_301$Metabolites %in% data_0.05$Metabolites,],starts_with("A1")), 
                         select(data_301[data_301$Metabolites %in% data_0.05$Metabolites,],starts_with("A2")))
RMETA2::savexlsx1(data_301heatmap,name = "heatmap.xlsx",sheet = "heatmap")
# heatmap
RMETA2::auto_alldraw(name="heatmap.xlsx",needgroup=NULL,drawwhat="heatmap",row=T,col=F,datasave = F,range = NULL)
#corr
RMETA2::auto_alldraw(name="heatmap.xlsx",drawwhat="correlation",datasave = T)


#####4
VOLCANO <- openxlsx::read.xlsx("volcano.xlsx",rowNames=T)
auto_volcano(VOLCANO,drawFC=2, name = "volcano")

#####5
##volcano
library(ggplot2)
vip_fc <- openxlsx::read.xlsx("差异代谢物(未筛选).xlsx")
vip_FC <- select(vip_fc, Metabolites,VIP, `log2(FC)`) %>% filter(VIP<10)
vip_FC$Class[vip_FC$VIP<=1] <- "not-significant"
vip_FC$Class[vip_FC$VIP>1&vip_FC$VIP<10&vip_FC$`log2(FC)`>=1] <- "up-regulated"
vip_FC$Class[vip_FC$VIP>1&vip_FC$VIP<10&abs(vip_FC$`log2(FC)`)<1] <- "VIP>1&|FC|<2"
vip_FC$Class[vip_FC$VIP>1&vip_FC$VIP<10&vip_FC$`log2(FC)`< -1]  <- "down-regulated"
p <- ggplot2::ggplot(data = vip_FC, mapping = aes(x = `log2(FC)`, y = VIP))     # 基础图层，不出现图形元素
p <- p + geom_point(aes(color = `Class`),size=2,alpha=0.8) +                                              
  scale_colour_manual(values=c("not-significant"="grey75","VIP>1&|FC|<2"="#CC9933", "down-regulated"="#66FF66", "up-regulated"= "#FF3300"))+
  labs(title="volcano", x="Log2(Fold Change)", y="VIP Score") +                              # 添加主题，x，y轴名称
  theme_light() +
  xlim(-9, 9)+
  theme(panel.grid =element_blank(), plot.title=element_text(hjust=0.5)
        ) +                # legend.position = "none"对主题进行修改，将中间网格消除，将主题名称移到正中间
  geom_hline(yintercept=1,linetype=2) +                                        # 添加水平、垂直分界线
  # geom_vline(xintercept=c(-1,1),linetype=2) +
  annotate("text", x = 7.5, y=0.75,label="VIP=1") 
ggsave(paste0("Volcano.pdf"), p, width = 9, height = 7)#, device = "pdf"
ggsave(filename = paste0("Volcano.tiff"), plot = p, width = 9, height = 7, dpi = 300)
statue = file.rename(paste0("Volcano.tiff"), paste0("Volcano.tif"))



#####6 
data_diff <- openxlsx::read.xlsx("差异代谢物.xlsx")
data_diff_heatmap <- cbind(select(data_diff,Metabolites), select(data_diff,starts_with("A1")), select(data_diff,starts_with("A2")))
RMETA2::savexlsx1(data_diff_heatmap, name = "heatmap.xlsx",sheet = "heatmap")
RMETA2::auto_alldraw(name="heatmap.xlsx",drawwhat="correlation",datasave = T)

# network
RMETA2::auto_alldraw(name="ID.xlsx",drawwhat = "network", p=NA, opacityNoHover = 1,species="hsa",fontSize = 14)
RMETA2::auto_alldraw(name="ID.xlsx",drawwhat = "network", p=NA, opacityNoHover = 0,species="hsa",fontSize = 14)

# 箱线图
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# ###################################################
# #导入数据
# #使用并列箱线图叠加散点图表示四缸，六缸，八缸发动机对每加仑汽油行驶的英里数
# #Group：为因子变量，表示汽车发动机的缸数
# #Value：为连续变量，表示对每加仑汽油行驶的英里数
# attach(mtcars)
# Data1 = data.frame(Group = mtcars$cyl,Value = mtcars$mpg)
# Data1$Group = factor(Data1$Group,levels = c(4,6,8))
# ###################################################

data_8 <- openxlsx::read.xlsx("经likelihood ratio test调整的混合效应识别差异代谢物数据矩阵—17种.xlsx")
data_box <- cbind(select(data_8,Metabolites), select(data_8,Class),
                  data.frame("low exposure" = log10(select(data_8,starts_with("A1")) %>% rowMeans()), check.names = F) , 
                  data.frame("high exposure" = log10(select(data_8,starts_with("A2")) %>% rowMeans()), check.names = F))
library(reshape2)
data_box_long <- reshape2::melt(data_box, id=c("Metabolites", "Class"), value.name="value") # Log Peak instensity
data_box_long <- transform(data_box_long, dist= pair(data_box_long$Class),scat_adj=ifelse(variable == "low exposure", -0.2, 0.2))

pair <- function(data=data, changedata=NULL){
  factordata <- factor(unique(data_box_long[,"Class"]))
  ord <- factordata[order(factordata)]
  if (is.null(changedata)){
    df_pair <- data.frame(ord, index=1:length(ord))
  }else{
  df_pair <- data.frame(ord, changedata)
  }
  change <- df_pair$index[match(data, df_pair$ord)]
  return(change)
}
pair(data_box_long$Class)
order(unique(data_box_long[,"Class"]))


#使用ggplot2包生成箱线图
p <- ggplot(data_box_long, aes(Class, value))+
    geom_boxplot(outlier.size = 0,
                 aes(color = factor(variable)),
                 position = position_dodge(1),
                 size=0.4)+
    geom_jitter(aes(dist+scat_adj, value, fill = factor(variable)),
                position = position_jitter(width=0.2, height = 0),
                shape=21,
                size=2)+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))+
 
    xlab("")+
    ylab("Log Peak instensity")+
    guides(fill=guide_legend(title = ""))+
       theme_bw()+ #背景变为白色
       theme(legend.position="righttop", #不需要图例
            panel.grid.major = element_blank(), #不显示网格线
            panel.grid.minor = element_blank())

# 采用下面方案，画箱线图
library(ggpubr)
  p <- ggpubr::ggboxplot(data_box_long, x="Class", y="value",
                    color="variable", add = "jitter", palette = "jama",
                    xlab = "",ylab = "Log Peak instensity")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1,size = 10),#
        legend.title=element_blank(),
        legend.justification=c(1.25,1.25),
        legend.position=c(1,1),
        panel.grid=element_blank(),
        plot.margin=unit(rep(3,4),'lines'),
        line = element_line(colour = "black"), 
        legend.margin = margin(t = 0, r = 6, b = 6, 
                               l = 2, unit = "pt"))+
  ylim(1, 1.2*max(data_box_long$value))
  ggsave(paste0("箱线图.pdf"), p, width =10, height = 8)#, device = "pdf"
  ggsave(filename = paste0("箱线图.tiff"), plot = p, width = 10, height = 8, dpi = 300)
  

#使用ggplot2包生成箱线图
P1 <- ggplot(data_box_long,aes(x=Class,y=value,fill=variable))+ #”fill=“设置填充颜色
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的性
  geom_jitter(aes(fill=variable),width =0.2,shape = 21,size=2.5)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  scale_fill_manual(values = c("black", "red"))+  #设置填充的颜色
  # scale_color_manual(values=c("black","black","black"))+ #设置散点图的圆圈的颜色为黑色
  ggtitle("Car Milleage Data")+ #设置总的标题
  theme_bw()+ #背景变为白色
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank())+
  ylab("Miles Per Gallon")+xlab("Number of Cylinders") #设置x轴和y轴的标题


## PCA loading 热图
# 提取数据
# 混合效应识别差异代谢物数据矩阵-8个

library(dplyr)
data_diff <- openxlsx::read.xlsx("经likelihood ratio test调整的混合效应识别差异代谢物数据矩阵—17种+FC.xlsx")
data <- cbind(select(data_diff,Metabolites), 
              select(data_diff,starts_with("A1")), 
              select(data_diff,starts_with("A2")))

rownames(data) <- data[,1]
data <- t(data[,-1]) %>% apply(., MARGIN = 2, as.numeric)

# PCA
# data_ <- data %>% prcomp(., scale = TRUE)
# data_pca <- data %>% princomp(.,cor = FALSE)
# 
# library(gmodels)
# data_ <- data %>%  fast.prcomp(.,scale. = T)
# summary(data_)
# summary(data_pca)
# predict(data_pca, data)[1:5,]
# screeplot(data_, type='lines')
# biplot(data_)

library(psych)
# 最佳组成分确定
print(fa.parallel(data,fa='pc',n.iter = 100,show.legend = FALSE))

# Parallel analysis suggests that the number of factors =  NA  and the number of components =  6
pca <- principal(data, nfactors = 5, rotate = "varimax", scores = TRUE)
# pca
# pca$scores
# PCA loading值

loadings <- data.frame(PC1 = pca$loadings[,1], PC2=pca$loadings[,2],
                       PC3 = pca$loadings[,3], PC4=pca$loadings[,4],
                       PC5 = pca$loadings[,5])

openxlsx::write.xlsx(loadings,file = "PC loadings.xlsx", sheetName  = "PC loadings", rowNames =T)
# loadings <-  data.table::setorderv(loadings,order=-1)
# 画图
# library(pheatmap)
# pheatmap(loadings,
#          scale = "none",                    #矩阵按行进行标准化
#          cluster_cols = F,                  #按列聚类
#          cluster_rows = T,                  #按行聚类
#          fontsize_row = 3.8,                #行字体大小             
#          fontsize_col = 5.5,                #列字体大小
#          cellwidth = 25, cellheigh = 6,     #设置单格大小
#          show_colnames = TRUE,              #是否显示行名称
#          color = colorRampPalette(c("blue", "white","red"))(100),  #设置渐变色
#          border_color = F,
#          angle_col = 0
#          # legend_breaks = c(-1, 0, 1),
#          # legend_labels = c(-1, 0, 1),
#          )  



simple_heatmap <- function(data=loadings, groupname = NULL, 
                           color=NULL,legendname=NULL,
                           column_names_rot=0,
                           cluster_rows=T,
                           show_row_dend=F,
                           ...){
  library(ComplexHeatmap)
  if (is.null(color)){
    color<-colorRampPalette(c("blue", "white","red"))(1000)
  }else{
    print("color argument must be availble like \"blue\" or \"#0000FF\" ")
    color<-color
    }
  if (is.null(legendname)){
    legendname <- "PC loadings"
  }else{
    legendname <- legendname
  }
  if (is.matrix(data)){
    data <- data
  }else{
    data <- as.matrix(data)
  }
  
  
  
  
  # heatmap body,width and height, By default, all heatmap components have fixed width or height
  if (ncol(data)<6){
    width <- ncol(data)
  }else if(ncol(data)<12){
    width <- 0.5*ncol(data)
  }else if(ncol(data)>=12){
    width <- 6
  }
  if (nrow(data)<=6){
    height <- 0.618*nrow(data)
  }else if(nrow(data)<12){
    height <- 10
  }else if(nrow(data)>=12){
    height <- 12
  } 
  

  # row col marker max 
  maxrowchar <- max(nchar(rownames(data)))
  maxcolchar <- max(nchar(colnames(data)))
  row_names_max_width <- max_text_width(rownames(data))
  column_names_max_height <- max_text_width(colnames(data))
  
  # legend width and height
  legend_height <- unit(height*0.4,"cm")
  legend_width  <- unit(height*0.618*0.3,"cm")
  
    # the complete heatmap including all heatmap components (excluding the legends) 
  heatmap_width <-  unit(width+height/nrow(data)*maxrowchar+height*0.618*0.3, "cm")
  heatmap_height <- unit(height+as.numeric(column_names_max_height)/10, "cm")
  
  # row col fontsize
  row_fontsize <- height/nrow(data)*72.27/2.54*0.8
  col_fontsize <- width/ncol(data)*72.27/2.54*0.5
  
   
  #  save figure
  # library(Cairo)
  # library(showtext)
  groupname <- ifelse(is.null(groupname),"",paste0(groupname,"-"))
  ## PDF                   
  # Cairo::CairoPDF(file = paste0(groupname,legendname,"-heatmap.pdf"),width = 12, height = 12)
  pdf(file = paste0(groupname,legendname,"-heatmap.pdf"),
      width =as.numeric(heatmap_width)/2, height = as.numeric(heatmap_height)/2)
  # showtext::showtext_begin()
  ComplexHeatmap::Heatmap(data,
        col = color,
        name = legendname,
        column_names_rot = column_names_rot,
        cluster_rows = cluster_rows,
        show_row_dend = show_row_dend,
        cluster_columns =F,
        width = width,
        height = height,
        row_names_max_width = row_names_max_width,
        column_names_max_height = column_names_max_height,
        heatmap_width = heatmap_width,
        heatmap_height = heatmap_height,
        row_names_gp = gpar(fontsize = row_fontsize),
        column_names_gp = gpar(fontsize = col_fontsize),
        row_names_centered = F,
        column_names_centered = T,
        heatmap_legend_param = list(legend_height = legend_height,
                                    legend_width = legend_width)
        )
  # showtext::showtext_end()
 # p2 <- rowAnnotation(text = row_anno_text(rownames(data)), width = max_text_width(rownames(data)))
 # draw(p1+p2)
  dev.off()
  # dev.new()
  
  ## tiff
  # Cairo::CairoTIFF(filename = paste0(groupname,legendname,"-heatmap.tiff"),width = 600, height = 600)
  tiff(filename = paste0(groupname,legendname,"-heatmap.tiff"),
       width = as.numeric(heatmap_width)*72.27/2, height = as.numeric(heatmap_height)*72.27/2, 
       units = "px")
  ComplexHeatmap::Heatmap(data,
                          col = color,
                          name = legendname,
                          column_names_rot = column_names_rot,
                          cluster_rows = cluster_rows,
                          show_row_dend = show_row_dend,
                          cluster_columns =F,
                          width = width,
                          height = height,
                          row_names_max_width = row_names_max_width,
                          column_names_max_height = column_names_max_height,
                          heatmap_width = heatmap_width,
                          heatmap_height = heatmap_height,
                          row_names_gp = gpar(fontsize = row_fontsize),
                          column_names_gp = gpar(fontsize = col_fontsize),
                          row_names_centered = F,
                          column_names_centered = T,
                          heatmap_legend_param = list(legend_height = legend_height,
                                                      legend_width = legend_width)
  )
  dev.off()
}

simple_heatmap(data=loadings, groupname = "17个差异代谢物")
# cxtCompMN <- cor(t(data), PCAloading[,1], use = "pairwise.complete.obs")




