load('tcga_rna_417sample.Rdata')
load(file ='E:/R/jiahui_doctor/lasso  RSF  批量单因素/risk_group.Rdata')
group=merge(risk_group,m,by='sample')
group=group[,1:2]
library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
row.names(m)=m[,1]
m=m[,-1]
m=as.data.frame(t(m))

#m<- tibble::rownames_to_column(m,"sample")#读数据


m=as.matrix(m)
m=log2(m+1)
m$group=group[,1]
m=as.data.frame(t(m))

n=substr(m$sample,1,12)
m$sample=n
m=merge(m,sample,by='sample')
row.names(m)=m[,1]
m=m[,-1]
m=as.data.frame(t(m))



data=m
#row.names(data)=data[,1]
#data=data[,-1]
#data=as.data.frame(t(data)) #行为gene列为样本
write.table(data,file='data_estimate.txt',sep="\t",
            row.names = T,col.names = T,quote = F)  #先存为txt
#data=read.table('data_estimate.txt',sep='\t',header = T,quote = '')

##关键函数，生成GCT文件
filterCommonGenes(input.f = "data_estimate.txt",output.f = "gene.gct",id="GeneSymbol")
#通过GCT文件评估肿瘤样本纯度
estimateScore("gene.gct","estimat_score.gct",platform = "affymetr")
score_table=read.table("estimat_score.gct",skip = 2,header = 1)

rownames(score_table)=score_table[,1]
score=t(score_table)
colnames(score)=score[1,]
score_matrix=score[-c(1,2),]

#######绘图
score_matrix=as.data.frame(score_matrix)
#score_matrix=score_matrix[-1,]

score_matrix$group=group[,2]
write.table(score_matrix,'estimete_result.csv',col.names = T,row.names = F,sep = ',')
score_matrix=read.csv('estimete_result.csv',header = T,sep=',')

library(ggplot2)
library(ggpubr)
#StromalScore
ggplot(score_matrix, aes(group, StromalScore),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 1600)
#ImmuneScore
ggplot(score_matrix, aes(group, ImmuneScore),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 2800)
#ESTIMATEScore
ggplot(score_matrix, aes(group, ESTIMATEScore),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 3200)
#TumorPurity
ggplot(score_matrix, aes(group, TumorPurity),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 1.05)  #加p值

