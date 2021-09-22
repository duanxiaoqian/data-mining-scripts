# LM2020-0005
setwd("Z:/项目数据/LCMS/2020/LM2020-0005/售后/售后1")
library(openxlsx)
library(dplyr)

d1_ve <- read.xlsx("差异代谢物.xlsx", sheet=1)
d3_ve <- read.xlsx("差异代谢物.xlsx", sheet=2)
d5_ve <- read.xlsx("差异代谢物.xlsx", sheet=3)
d3_d1 <- read.xlsx("差异代谢物.xlsx", sheet=4)
d5_d1 <-  read.xlsx("差异代谢物.xlsx", sheet=5)
d5_d3 <-  read.xlsx("差异代谢物.xlsx", sheet=6)

total <-  read.xlsx("数据矩阵.xlsx", sheet=1)
# 并集

diff_union <- unique(c(d1_ve$Metabolites, d3_ve$Metabolites, d5_ve$Metabolites,
                     d3_d1$Metabolites, d5_d1$Metabolites, d5_d3$Metabolites))
diff_union_data <- total[total$Metabolites %in% diff_union,] %>% select(., 19:42)
diff_union_data_final <- data.frame(Metabolites=diff_union,
                                    Vehicle=select(diff_union_data, starts_with("Vehicle")) %>% rowMeans(),
                                    Day1=select(diff_union_data, starts_with("Day1")) %>% rowMeans(),
                                    Day3=select(diff_union_data, starts_with("Day3")) %>% rowMeans(),
                                    Day5=select(diff_union_data, starts_with("Day5")) %>% rowMeans())
RMETA2::savexlsx1(diff_union_data_final, "diff_STR.xlsx", sheet="diff_union")
# intersect
diff_intersect <- multi_intersect(List = list(d1_ve$Metabolites, d3_ve$Metabolites, d5_ve$Metabolites,
                            d3_d1$Metabolites, d5_d1$Metabolites, d5_d3$Metabolites))
diff_intersect_data <- total[total$Metabolites %in% diff_intersect,]
diff_intersect_data_final <- data.frame(Metabolites=diff_intersect,
                                        Vehicle=select(diff_intersect_data, starts_with("Vehicle")) %>% rowMeans(),
                                        Day1=select(diff_intersect_data, starts_with("Day1")) %>% rowMeans(),
                                        Day3=select(diff_intersect_data, starts_with("Day3")) %>% rowMeans(),
                                        Day5=select(diff_intersect_data, starts_with("Day5")) %>% rowMeans())

RMETA2::savexlsx1(diff_intersect_data_final, "diff_STR.xlsx", sheet="diff_intersect")
RMETA2::dSTEM(filename = "diff_STR.xlsx",sheet = 1,cluster = 15)
dSTEM(filename="diff_STR.xlsx", sheet=2,
      circle_num=2, max_cluster=50, cluster_num=NULL,
      filter_std=FALSE, standardise=T, mode="metabolites",
      color=NULL
      )
# RMETA2::dSTEM(filename = "diff_STR.xlsx",sheet =2)




setwd("Z:/项目数据/LCMS/2020/LM2020-0005/售后/售后1")

