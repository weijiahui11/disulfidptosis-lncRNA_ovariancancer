library(maftools)
dir.path="E:/1.组学数据/OV/TCGA_OV_maf"
all.maf=list.files(path = dir.path,pattern = ".gz",
                   full.names = T,recursive = T)
all.maf[1:3]
maf.list=lapply(all.maf,data.table::fread,
                sep='\t',
                header=T,
                skip=7)
maf=do.call(rbind,maf.list)
dim(maf)
maf1=maf  ##备用
maf=read.maf(maf1) 

tmb = function(maf, captureSize = 50, logScale = TRUE){
  
  if(!is(object = maf, class2 = "MAF")){
    stop("Input must be an MAF object")
  }
  
  maf.mutload = getSampleSummary(maf)[,.(Tumor_Sample_Barcode, total)]
  maf.mutload[,total_perMB := total/captureSize]
  maf.mutload[,total_perMB_log := log10(total_perMB)]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  
  medload = median(maf.mutload[,total_perMB], na.rm = TRUE)
  
  par(mar = c(2, 4, 2, 0))
  pointcol = "#009688"
  medlinecol = "#FF5722"
  if(logScale){
    yat = pretty(range(maf.mutload[total != 0][,total_perMB_log]))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
    abline(h = log10(medload), lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB_log, pch = 19, col = pointcol)
    title (main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2,at = yat, las = 2)
    mtext(text = "TMB/MB (log10)", side = 2, line = 2.5)
  }else{
    yat = pretty(range(maf.mutload$total_perMB))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
    abline(h = medload, lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB, pch = 19, col = pointcol)
    title (main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2,at = yat, las = 2)
    mtext(text = "TMB/MB", side = 2, line = 2.5)
  }
  
  
  maf.mutload
}
TMB = tmb(maf = maf)



a=maf@data
sample=data.frame(sam=a$Tumor_Sample_Barcode)
names(sample)='Tumor_Sample_Barcode'
sample$sample=substr(sample$Tumor_Sample_Barcode,1,12)
sample=merge(sample,risk_group,by='sample')
##分高低risk组
high=sample[which(sample$group=='high'),]
high=high[!duplicated(high$Tumor_Sample_Barcode),]
high=as.data.frame(high[,2])
names(high)='Tumor_Sample_Barcode'

low=sample[which(sample$group=='low'),]
low=low[!duplicated(low$Tumor_Sample_Barcode),]
low=as.data.frame(low[,2])
names(low)='Tumor_Sample_Barcode'


####绘图
h=a
l=a
h=merge(a,high,by='Tumor_Sample_Barcode')

##加入临床分组信息
group=maf@clinical.data
group_h=merge(high,surv,by='Tumor_Sample_Barcode')
group_h=merge(group,group_h,by='Tumor_Sample_Barcode') ##必须将group放在前面 不然格式不对
group_l=merge(low,surv,by='Tumor_Sample_Barcode')
group_l=merge(group,group_l,by='Tumor_Sample_Barcode')

library(maftools)
maf@data=h
maf@clinical.data=group_h
plotmafSummary(maf=maf,rmOutlier = T,addStat = 'median',dashboard = T)
#改颜色
fabcolors = "#D53E4F"
names(fabcolors) = c('high')
fabcolors = list(group = fabcolors)
oncoplot(maf =maf, fontSize = 0.6 ,showTumorSampleBarcodes = F ,clinicalFeatures = 'group',titleFontSize=1.2,legendFontSize=0.7,removeNonMutated=F,writeMatrix=T,
         annotationFontSize=0.7,sortByAnnotation = TRUE,
         annotationColor = fabcolors)
#oncoplot(maf = maf)

l=merge(a,low,by='Tumor_Sample_Barcode')
maf@data=l
maf@clinical.data=group_l
plotmafSummary(maf=maf,rmOutlier = T,addStat = 'median',dashboard = T)
#改颜色
fabcolors = "#3288BD"
names(fabcolors) = c('low')
fabcolors = list(group = fabcolors)
oncoplot(maf =maf, fontSize = 0.6 ,showTumorSampleBarcodes = F ,clinicalFeatures = 'group',titleFontSize=1.2,legendFontSize=0.7,removeNonMutated=F,writeMatrix=T,
         annotationFontSize=0.7,sortByAnnotation = TRUE,
         annotationColor = fabcolors)
#oncoplot(maf = maf)

###加入临床分组信息
group_h=merge(high,surv,by='Tumor_Sample_Barcode')


###相关性图
par(oma=c(3,4,5,1))
par(oma = c(3, 4, 5, 1))
somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.1),fontSize = 0.7)

###
maf@data[["TumorVAF"]] <- maf@data$t_alt_count / maf@data$t_depth
plotVaf(maf = maf, vafCol = 'TumorVAF')



##计算tmb
tmb = function(maf, captureSize = 50, logScale = TRUE){
  
  if(!is(object = maf, class2 = "MAF")){
    stop("Input must be an MAF object")
  }
  
  maf.mutload = getSampleSummary(maf)[,.(Tumor_Sample_Barcode, total)]
  maf.mutload[,total_perMB := total/captureSize]
  maf.mutload[,total_perMB_log := log10(total_perMB)]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]
  
  medload = median(maf.mutload[,total_perMB], na.rm = TRUE)
  
  par(mar = c(2, 4, 2, 0))
  pointcol = "#009688"
  medlinecol = "#FF5722"
  if(logScale){
    yat = pretty(range(maf.mutload[total != 0][,total_perMB_log]))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
    abline(h = log10(medload), lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB_log, pch = 19, col = pointcol)
    title (main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2,at = yat, las = 2)
    mtext(text = "TMB/MB (log10)", side = 2, line = 2.5)
  }else{
    yat = pretty(range(maf.mutload$total_perMB))
    plot(NA, xlim = c(0, nrow(maf.mutload)), ylim = range(yat), axes = FALSE, xlab = NA, ylab = NA)
    abline(h = yat, lty = 2, lwd = 1, col = 'gray90')
    abline(h = medload, lty = 2, lwd = 1, col = medlinecol)
    points(x = 1:nrow(maf.mutload), y = maf.mutload$total_perMB, pch = 19, col = pointcol)
    title (main = "Mutation Burden", sub = paste0("Median: ", medload, "/MB"), line = 0, adj = 0)
    axis(side = 2,at = yat, las = 2)
    mtext(text = "TMB/MB", side = 2, line = 2.5)
  }
  
  
  maf.mutload
}
TMB = tmb(maf = maf)
t=TMB[,1:2]
#t=merge(t,sample,by='Tumor_Sample_Barcode')
s=sample
s=s[!duplicated(s$Tumor_Sample_Barcode),]
t=merge(t,s,by='Tumor_Sample_Barcode')
#TMB分高低组
n=as.matrix(t[,2])
median(n[,1])
t[,2]<-ifelse(t[,2]<=57,"L","H")
##加载生存数据
load('E:/R/jiahui_doctor/lasso  RSF  批量单因素/surv.RData')
surv$sample=row.names(surv)
t=merge(t,surv,by='sample')
write.table(t,'tmb.csv',col.names = T,row.names = F,sep = ',')
t=read.csv('tmb.csv',header = T,sep=',')
surv=t  ##不能用t绘制生存曲线 会报错

library(ggpubr)
library(ggplot2)
library(rms)  # 画列线图
library(survival) # 生存分析包
library(survminer)
f<-survfit(Surv(time=surv$time,event = surv$os)~TMB,
           data=surv)
ggsurvplot(f,data = surv,
           #title  = "PI值高低组患者的生存曲线",  
           xlab = "Time(days)", 
           ylab = "OS",
           font.main = c(12, "bold", "darkblue"), #字体大小，样式和颜色
           font.x = c(10, "darkred"), #x轴标签的字体大小，斜体和颜色
           font.y = c(10, "darkred"),#y轴的字体大小，颜色
           font.tickslab = c(10, "plain", "darkgreen"), #更改刻度标签字体大小，样式颜色
           ggtheme = theme_bw(),
           pval = TRUE,#增加p值
           pval.coord = c(800, 0.1),#p值位置坐标
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


###小提琴图
library(ggplot2)
library(ggpubr)
names(t)[2]='TMB'
ggplot(t, aes(group, TMB),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 300)