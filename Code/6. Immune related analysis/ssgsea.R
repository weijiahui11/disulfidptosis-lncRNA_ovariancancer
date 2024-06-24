library(ggplot2)
library(tinyarray)
library(dplyr)
library(Hmisc)
library(pheatmap)
#BiocManager::install("GSVA")
library(GSVA)

geneset = rio::import("E:/魏佳慧/宫颈癌/cellMarker.xlsx")
geneset = split(geneset$Metagene,geneset$`Cell type`) #加载后转化为list
lapply(geneset[1:3], head) 

load('E:/R/jiahui_doctor/卵巢癌处理/tcga_rna.RData')
mrna_count= tibble::rownames_to_column(m,"sample")
mrna_count=merge(mrna_count,sample,by='sample')
row.names(mrna_count)=mrna_count[,1]
mrna_count=mrna_count[,-1]
mrna_count=as.data.frame(t(mrna_count))
mrna_count=log2(mrna_count+1)

#数据行为gene列为样本 并转为矩阵
ssgsea=mrna_count
#row.names(ssgsea)=ssgsea[,1]
#ssgsea=ssgsea[,-1]
#ssgsea=as.data.frame(t(ssgsea))
ssgsea<- as.matrix(ssgsea)
#进行gsva分析
re <- gsva(ssgsea, geneset, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE) 
#group=load(file='E:/R/jiahui_doctor/批量单因素/group.Rdata')
gp=group
gp=gp[,2]
draw_boxplot(re,gp,color=c("#1d4a9b","#e5171a"))
