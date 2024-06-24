library(pRRophetic)
#data(bortezomibData)
#dim(exprDataBortezomib)


load('E:/R/jiahui_doctor/卵巢癌处理/mrna_fpkm.RData')
load("E:/R/jiahui_doctor/lasso  RSF  批量单因素/risk_group.RData")
mrna_count=merge(risk_group,mrna_count,by='sample')
risk_group=mrna_count[,1:2]
mrna_count=mrna_count[,-2]
row.names(mrna_count)=mrna_count[,1]
mrna_count=mrna_count[,-1]
mrna_count=as.data.frame(t(mrna_count))  ##行gene  列样本
mrna_count=log2(mrna_count+1)
mrna_count=as.matrix(mrna_count)


predictedPtype <- pRRopheticPredict(testMatrix=mrna_count, 
                                    drug="EHT 1864",  #不同药物
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2016")  #2016药物种类多
head(predictedPtype) ##结果
predictedPtype=as.data.frame(predictedPtype)
#row.names(risk_group)=risk_group[,1]
predictedPtype=cbind(predictedPtype,risk_group)
names(predictedPtype)[1]='value'

library(ggplot2)
library(ggpubr)

ggboxplot(predictedPtype, x="group",y="value",fill="group",alpha=0.3,palette = "lancet",
          ylab="EHT 1864 Sensitivity (IC50)",
          xlab="Risk Group")+
  stat_compare_means(method = "t.test",label.y =5)+  ##修改p位置
  theme(legend.position = "none")




