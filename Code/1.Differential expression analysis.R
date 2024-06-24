load('E:/1.组学数据/gtex联合tcga/tcga_lnc.RData')
load('E:/1.组学数据/gtex联合tcga/gtex_lnc.RData')
#tcga去重
row.names(tcga_lnc)=tcga_lnc[,1]
tcga_lnc=tcga_lnc[,-1]
tcga_lnc=as.data.frame(t(tcga_lnc))
tcga_lnc<- tibble::rownames_to_column(tcga_lnc,"sample")
a=substr(tcga_lnc$sample,1,12)
tcga_lnc$sample=a
tcga_lnc=tcga_lnc[which(!duplicated(tcga_lnc$sample)),]

data=merge(tcga_lnc,gtex_lnc,by='gene')
row.names(data)=data[,1]
data=data[,-1]

library(dplyr)
library('gplots')
library('limma')

group=as.data.frame(colnames(data))
names(group)='sample'
group$group=group$sample
write.table(group,'group.csv',col.names = T,row.names = F,sep=',')
group=read.csv('group.csv',header=T,sep=',')

foldChange=1 #fold change=1意思是差异是两倍
padj=0.05#p

group <- group[,2]#表示分组那一列
design <- model.matrix(~0+factor(group))
colnames(design) <- c("normal","tumor")  #要改名 不然会报错
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
diffSig = nrDEG[(nrDEG$adj.P.Val < padj & 
                   (nrDEG$logFC>foldChange | nrDEG$logFC<(-foldChange))),]

up=nrDEG[(nrDEG$adj.P.Val < padj&nrDEG$logFC>foldChange),]
down=nrDEG[(nrDEG$adj.P.Val < padj&nrDEG$logFC<(-foldChange)),]
lnc_dif=as.data.frame(rownames(diffSig))
names(lnc_dif)='gene'
lnc_dif=merge(lnc_dif,tcga_lnc,by='gene')
save(lnc_dif,file = 'lnc_dif.RData')


###火山图 
library(EnhancedVolcano)
library(airway)
EnhancedVolcano(nrDEG,
                lab=rownames(nrDEG),
                x='logFC',
                y='adj.P.Val',
                xlim=c(-8,8),
                title='DESeq2_DEG',
                pCutoff = 5e-2,
                FCcutoff=1,
                pointSize = 1,
                labSize = 3.0,
                col=c('black','green','blue','red3'),
                shape = 8,
                colAlpha = 1,
                legend = c('NS','LogFC','P value','P value & LogFC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                axisLabSize = 10)









######热图
library(pheatmap)
library(RColorBrewer)
group=read.csv('group.csv',header=T,sep=',')
d=as.data.frame(group$group)
names(d)='group'
color<-colorRampPalette(c('#436eee','white','#EE0000'))(100)
row.names(d)=colnames(data)


#data1=data
#data1=log2(data1+0.01)

#pheatmap(data1,show_rownames = F,show_colnames = F,
             #cluster_rows = F,cluster_cols = F,
             #annotation_col = d,color = color,
             #legend = T ,
             #scale = "row",border_color = black
#)
