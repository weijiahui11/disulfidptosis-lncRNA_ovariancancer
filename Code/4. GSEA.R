load(file = "E:/R/jiahui_doctor/卵巢癌处理/tcga_rna.Rdata")
load(file = 'E:/R/jiahui_doctor/lasso  RSF  批量单因素/risk_group.Rdata')

row.names(m)=m[,1]
m=m[,-1]
m=as.data.frame(t(m))
m<- tibble::rownames_to_column(m,"sample")
a=substr(m$sample,1,12)
m$sample=a
sample=as.data.frame(risk_group[,1])
names(sample)='sample'
m=merge(m,sample,by='sample')
#sample=as.data.frame(m[,1])
group=merge(risk_group,sample,by='sample')
row.names(m)=m[,1]
m=m[,-1]
m=as.data.frame(t(m))
m=log2(m+1)


#######差异分析
library(dplyr)
library('gplots')
library('limma')

foldChange=0.5#fold change=1意思是差异是两倍
padj=0.05#
group <-group[,2]#表示分组那一列
design <- model.matrix(~0+factor(group))
colnames(design) <- c("high","low")

fit <- lmFit(m,design)
cont.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)  


##要用全基因做GSEA分析 不能只选择差异基因
#######GSEA
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

df=data.frame(SYMBOL=rownames(nrDEG),
              logFC=nrDEG$logFC)
rownames(df)=rownames(nrDEG)
df_id<-bitr(df$SYMBOL, #转换的列是df数据框中的SYMBOL列  unique(df$SYMBOL)
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型  c("ENTREZID")
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
df_all<-merge(df,df_id,by="SYMBOL",all=F)#使用merge合并
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$logFC #把foldchange按照从大到小提取出来
names(gene_fc) <- df_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(gene_fc)

library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

k=gseKEGG(geneList = gene_fc,
          organism = 'hsa',
          nPerm=10000,
          pvalueCutoff = 0.05,
          pAdjustMethod = "none")#p值校正方法)
class(k)          
colnames(k@result)
kegg_result=as.data.frame(k)
rownames(k@result)[head(order(k@result$enrichmentScore))]
af=as.data.frame(k@result)
write.table(af,file = paste0("2.","all_GSEA.xls"),sep="\t",quote = F,col.names = T)

gseaplot2(k, base_size = 10,
          geneSetID = row.names(t)[c(1:8)],
          rel_heights = c(0.5, 0.11, 0.11),
          ES_geom = "line")