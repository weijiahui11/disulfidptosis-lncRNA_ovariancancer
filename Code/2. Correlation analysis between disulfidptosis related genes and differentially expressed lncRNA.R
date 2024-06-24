#行为样本 列为gene
load('E:/R/jiahui_doctor/TCGA_GTEX联合2/lnc_dif.RData')
load('E:/1.组学数据/gtex联合tcga/tcga_liu.RData')
row.names(lnc_dif)=lnc_dif[,1]
lnc_dif=lnc_dif[,-1]
lnc_dif=as.data.frame(t(lnc_dif))
row.names(tcga_liu)=tcga_liu[,1]
tcga_liu=tcga_liu[,-1]
tcga_liu=as.data.frame(t(tcga_liu))

library(psych)
cor=corr.test(lnc_dif,tcga_liu)
cor_r=cor$r
cor_p=cor$p
library(reshape2)
cor_r=melt(as.matrix(cor_r))
cor_p=melt(as.matrix(cor_p))
melt_data=merge(cor_r,cor_p,by=c("Var1","Var2"),sort = F)
colnames(melt_data)=c("Var1","Var2","r","p")
lnc_m_cor = melt_data[(melt_data$p < 0.001 &(melt_data$r>0.3 | melt_data$r<(-0.3))),]
lnc_liu=as.data.frame(lnc_m_cor[,1])
names(lnc_liu)='lnc'
lnc_liu$gene=lnc_liu$lnc
lnc_liu=lnc_liu[which(!duplicated(lnc_liu$lnc)),]
lnc_liu=as.data.frame(lnc_liu[,-1])
names(lnc_liu)='gene'

###重新跑 不需要跑
lnc_dif=as.data.frame(t(lnc_dif))
lnc_dif<- tibble::rownames_to_column(lnc_dif,"gene")
#####

lnc_liu_exp=merge(lnc_liu,lnc_dif,by='gene')  ###高相关性lnc
row.names(lnc_liu_exp)=lnc_liu_exp[,1]
lnc_liu_exp=lnc_liu_exp[,-1]
#save(lnc_liu_exp,file = 'lnc_liu_exp.RData')
###加载生存数据  确定样本量
load('E:/R/jiahui_doctor/lasso  RSF  批量单因素/surv.RData')
lnc_liu_exp=as.data.frame(t(lnc_liu_exp))
lnc_liu_exp<- tibble::rownames_to_column(lnc_liu_exp,"sample")
a=substr(lnc_liu_exp$sample,1,12)
lnc_liu_exp$sample=a
lnc_liu_exp=merge(lnc_liu_exp,sample,by='sample')
row.names(lnc_liu_exp)=lnc_liu_exp[,1]
lnc_liu_exp=lnc_liu_exp[,-1]
save(lnc_liu_exp,file = 'lnc_liu_exp.RData')
