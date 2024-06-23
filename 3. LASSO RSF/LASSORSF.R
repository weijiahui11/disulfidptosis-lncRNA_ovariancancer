
###########读取生存、临床数据
#surv=read.table('E:/1.组学数据/OV/ov_ucsc_surv.tsv',header=T,sep='\t')
#surv=surv[which(!duplicated(surv$sample)),]
load('surv.RData')
load('E:/R/jiahui_doctor/lnc_liu_exp.RData')
data=cbind(lnc_liu_exp,surv)


####不需要跑
lnc_liu_exp=as.data.frame(t(lnc_liu_exp))
lnc_liu_exp<- tibble::rownames_to_column(lnc_liu_exp,"sample")
a=substr(lnc_liu_exp$sample,1,12)
lnc_liu_exp$sample=a
lnc_liu_exp=lnc_liu_exp[which(!duplicated(lnc_liu_exp$sample)),]
sample=as.data.frame(surv[,1])
names(sample)='sample'
lnc_liu_exp=merge(lnc_liu_exp,sample,by='sample')
sample=as.data.frame(lnc_liu_exp[,1])
names(sample)='sample'
surv=merge(surv,sample,by='sample')
#########


data=merge(lnc_liu_exp,surv,by='sample')
#row.names(lnc_liu_exp_surv)=lnc_liu_exp_surv[,1]
#lnc_liu_exp_surv=lnc_liu_exp_surv[,-1]
#data=lnc_liu_exp_surv



#data=cbind(lnc_liu_exp,surv)
###lasso cox
library(glmnet) ##Lasso回归
library(survival)
set.seed()
x = as.matrix(data[, 1:2467])
y = data.matrix(Surv(data$time,data$os))
fitcv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=5)
plot(fitcv)
coef(fitcv, s="lambda.min")

fit<- glmnet(x,y,family="cox", alpha=1,nfolds=5)
plot(fit)

c=as.data.frame(rownames(coef(fitcv, s = 'lambda.min'))[coef(fitcv, s = 'lambda.min')[,1]!= 0])
names(c)='gene'

#lnc_liu_exp=as.data.frame(t(lnc_liu_exp))
#lnc_liu_exp<- tibble::rownames_to_column(lnc_liu_exp,"gene")




c=merge(lnc_liu_exp,c,by="gene")
row.names(c)=c[,1]
c=c[,-1]
c=as.data.frame(t(c))
c$time=surv$time
c$status=surv$os


###随机生存森林
library(randomForestSRC)
library(survcomp)
RSF <- rfsrc(Surv(time, status) ~ ., data=c,importance =TRUE)
#importance <- RSF[["importance"]]
plot(RSF)
risk=predict(RSF,newdata = c)
chf <- rowSums(risk$chf)
cindex <- concordance.index(chf,surv.time = c$time
                            , surv.event = c$status)
cindex1 <- cindex$c.index


###多因素cox
l=c("CTD-2371O3.2","CTD-3193O13.11","RP11-872J21.3","AC074286.1",
    "CTB-171A8.1","IGFBP7-AS1","TLR8-AS1","CH507-254M2.2","RP11-102K13.5",
    "RP11-670E13.6")
duo=c[,which(colnames(c)%in%l)]
duo$time=c$time
duo$status=c$status


library(survival)
cox=coxph(Surv(time,status)~.,data=duo)
summary(cox)
res <- predict(cox,duo)
sum.surv<-summary(cox)
c_index<-sum.surv$concordance
c_index


### K-M生存曲线  ROC
#row.names(surv)=surv[,1]
#surv=surv[,-1]


#surv$risk=res
#risk_group=data.frame(sample=row.names(surv),group=surv$risk)
#save(risk_group,file = 'risk_group.RData')


##########用截断值绘制生存曲线
surv$risk=chf
#计算最佳截点
res.cut <- surv_cutpoint(surv, #数据集
                         time = "time", #生存状态
                         event = "os", #生存时间
                         variables = c("risk") )#需要计算的数据列名
summary(res.cut)#查看最佳截断值
plot(res.cut, "risk", palette = "npg")
#根据截点分类数据
res.cat <- surv_categorize(res.cut)
head(res.cat)

f<-survfit(Surv(time=res.cat$time,event = res.cat$os)~risk,
           data=res.cat)
ggsurvplot(f,data = res.cat,
           #title  = "PI值高低组患者的生存曲线",  
           xlab = "Time(days)", 
           ylab = "OS",
           font.main = c(12, "bold", "darkblue"), #字体大小，样式和颜色
           font.x = c(10, "black"), #x轴标签的字体大小，斜体和颜色
           font.y = c(10, "black"),#y轴的字体大小，颜色
           font.tickslab = c(10, "plain", "darkgreen"), #更改刻度标签字体大小，样式颜色
           ggtheme = theme_bw(),
           pval = TRUE,#增加p值
           pval.coord = c(1000, 0.1),#p值位置坐标
           #pval.coord = c(500, 0.1),#p值位置坐标
           pval.size =4,#p值字体大小
           conf.int = F, #置信区间
           palette = "lancet", # 自定义调色板
           legend.title = "risk",
           legend = c(0.12, 0.3),
           pval.method=TRUE,  #增加long rank检验
           pval.method.size=4, 
           pval.method.coord=c(0,0.1) , # 指定图例位置
           size = 1,  #更改线号
           linetype = "strata" #按组更改线型（即“分层”）直线与虚线分为两组
)

risk_group=data.frame(sample=row.names(res.cat),group=res.cat$risk)
save(risk_group,file = 'risk_group.RData')



###时间依赖性roc

surv$risk=chf
library(timeROC)
time_roc_res <- timeROC(
  T = surv$time,
  delta = surv$os,
  marker = surv$risk,
  cause = 1,
  weighting="marginal",
  times = c(1 * 365, 3 * 365, 5 * 365),
  ROC = TRUE,
  iid = TRUE
)
#计算AUC值及其置信区间
time_roc_res$AUC
#查看AUC的95%置信区间
confint(time_roc_res, level = 0.95)$CI_AUC
#绘制time-dependent ROC曲线
plot(time_roc_res, time=1 * 365, col = "red", title = FALSE)  
plot(time_roc_res, time=3 * 365, add=TRUE, col="blue") 
plot(time_roc_res, time=5 * 365, add=TRUE, col="green") 
legend("bottomright",c("1 Years" ,"3 Years", "5 Years"),
       col=c("red", "blue", "green"), lty=1, lwd=2)

#也可以通过修改在再美观点，
time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

#通过plotAUCcurve函数绘制不同时间节点的AUC曲线及其置信区间，也可将多个ROC曲线的AUC值放在一起绘制
plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
legend("bottomright", "mayoscore4", col = "red", lty=1, lwd=2)

#对于上述timeROC的结果，如3年ROC曲线的约登指数
surv$risk[which.max(time_ROC_df$TP_3year - time_ROC_df$FP_3year)]
