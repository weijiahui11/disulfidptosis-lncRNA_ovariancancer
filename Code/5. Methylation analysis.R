#mrna=as.data.frame(t(mrna))
load('E:/1.组学数据/gtex联合tcga/tcga_mRNA.RData')
#names(mrna)[1]='gene'
G=read.csv('甲基化m1a.csv',header = T,sep=',')
G=merge(G,m,by='gene')
row.names(G)=G[,1]
G=G[,-1]
G=as.data.frame(t(G))
load("E:/R/jiahui_doctor/lasso  RSF  批量单因素/lnc_8.RData")

G<- tibble::rownames_to_column(G,"sample")
n=substr(G$sample,1,12)
G$sample=n
sample=as.data.frame(row.names(lnc_8))
names(sample)='sample'
G=merge(sample,G,by='sample')

#n<-gsub("\\.","-",G$sample)
#G$sample<-n
#n=read.csv('a.csv',header = T,sep=',')
#G=merge(G,n,by='sample')
#n=substr(G$sample,1,12)
#G$sample=n
#G=G[which(!duplicated(G$sample)),]
#G=merge(G,c,by='sample')
row.names(G)=G[,1]
G=G[,-1]

#G=as.matrix(G)
#write.table(G,'G.csv',row.names = F,col.names = T,sep=',')
#G=read.csv('G.csv',header = T,sep=',')
#l=read.csv('12个lncRNA.csv',header = T,sep=',')
#lncRNA_diff=as.data.frame(t(lncRNA_diff))
#lncRNA_diff<- tibble::rownames_to_column(lncRNA_diff,"gene")
#l=merge(l,lncRNA_diff,by='gene')
#row.names(l)=l[,1]
#l=l[,-1]
#l=as.data.frame(t(l))
#l=as.matrix(l)

library(psych)
R=corr.test(lnc_8,G)
cor_r=R$r
cor_p=R$p

library(reshape2)
cor_r=melt(as.matrix(cor_r))
cor_p=melt(as.matrix(cor_p))
melt_data=merge(cor_r,cor_p,by=c("Var1","Var2"),sort = F)
colnames(melt_data)=c("lnc","m","r","p")
#修改P值分类
write.table(melt_data,'甲基化相关性.csv',row.names = F,col.names = T,sep=',')
melt_data=read.csv('甲基化相关性.csv',header = T,sep=',') 

library(ggplot2)
p=ggplot(melt_data)+
  geom_point(aes(lnc,m,size=abs(r),fill=p),color='#999999',
             shape=21)+
  scale_fill_manual(values = c('burlywood','darkcyan','purple'))+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1))+
  xlab("")+
  ylab("")+
  guides(size=guide_legend(title = "Spearman's |r|",order = 1),
         fill=guide_legend(title="p value",order = 2),
         color=guide_legend(override.aes = list(size=10)))
p+ggtitle("m1A") +
  theme(plot.title = element_text(hjust=0.5))

