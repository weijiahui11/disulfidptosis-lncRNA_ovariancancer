library(CIBERSORT)
library(ggplot2)#绘图
library(pheatmap)#绘制热图
library(ggpubr)#堆积比例图
library(reshape2)#数据处理
library(tidyverse)#数据处理

# 准备数据:转录组+免疫细胞类型列表+分组数据 -------------------------------------
#行gene列样本
load("E:/R/jiahui_doctor/卵巢癌处理/mrna_fpkm.RData")
load("E:/R/jiahui_doctor/group.RData")
mrna_fpkm=merge(group,mrna_fpkm,by='sample')
mrna_fpkm=mrna_fpkm[,-2]

row.names(mrna_fpkm)=mrna_fpkm[,1]
mrna_fpkm=mrna_fpkm[,-1]
mrna_fpkm=as.data.frame(t(mrna_fpkm))
DEG_expr <- mrna_fpkm
DEG_expr=as.matrix(DEG_expr)
#group <- read.csv("group.csv",row.names = 1)


#group$samples=rownames(group)

##行gene列样本
#数据检验
boxplot(DEG_expr,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(DEG_expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)#判断是否进行log的标准
LogC
if (LogC) { 
  DEG_expr[which(DEG_expr <= 0)] <- NaN
  DEG_expr <- log2(DEG_expr) 
  print("log2 transform is finished")
} else {
  print("log2 transform is not needed")
}   

# 均数不一致需标准化 ---------------------------------------------------------------
library(limma)
DEG_expr=normalizeBetweenArrays(DEG_expr)
boxplot(DEG_expr,outline=FALSE, notch=F , las=2)

# 免疫细胞类型列表 ----------------------------------------------------------------
LM22_local <- read.table("LM22.txt",header = T,row.names = 1,sep = "\t")
#data(LM22)#调用cibersort包中自带的LM22数据，22种免疫细胞参考marker基因表达情况
#all(LM22==LM22_local)#判断官网下载的文件是否与cibersort包中自带的LM22数据一致

# 正式开始Cibersort免疫细胞浸润分析 ---------------------------------------------------
#运行CIBERSORT以估计样本中免疫细胞类型的比例
result <- cibersort(sig_matrix = LM22, mixture_file = DEG_expr, perm = 1, QN = TRUE)
#参数包括 sig_matrix、mixture_file、perm 和 QN。sig_matrix 参数是 CIBERSORT 软件包的内置数据，mixture_file 参数是待测样本的基因表达数据。perm 参数表示是否使用随机排列法，QN 参数是 TRUE 或 FALSE，表示是否使用质量归一化(芯片数据设置为T，测序数据就设置为F）。cibersort 函数会返回一个数据框，包含了各种免疫细胞类型的比例和基因表达数据的质量归一化结果。

result <- as.data.frame(result)
colnames(result)
write.csv(result,"cibersort_result.csv")
#`P-value`越小越可信；
#Correlation原表达矩阵乘以细胞占比后的数据矩阵与原表达矩阵的相关性
#RMSE 均方根误差，越小效果越好

# 结果可视化 -------------------------------------------------------------------
result1 <- result[,1:ncol(LM22)]
result1 <- result1[,apply(result1, 2, function(x){sum(x)>0})]#删除全是0的列
#在矩阵的列上应用函数：apply(X, 2, fun)，其中X是矩阵，fun是要应用的函数，2表示按列应用函数。
#热图
pdf(file="Heatmap.pdf",width = 10,height = 8)
pheatmap(result1,
         color = colorRampPalette(c("#4CB43C", "#FEFCFB", "#ED5467"))(100),
         border="skyblue",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = T,#是否展示行名
         show_colnames = T,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'row',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         cluster_col=T,#分别指定是否按列和行聚类。
         cluster_row=F,
         angle_col = "45",#指定列标签的角度。
         legend = F,#指定是否显示图例。
         legend_breaks=c(-3,0,3),#指定图例中显示的数据范围为-3到3。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)
dev.off()


# 盛夏的果实函数科普一： --------------------------------------------------------------
####在R语言中，melt函数是用来将一个数据集转化为“长格式”的。
###所谓的“长格式”是指将一个数据集的多个变量的值放在一个列中，
###而不是将每个变量放在不同的列中。这样就可以使用一列来记录变量的名称，
###另一列来记录变量的值。

# ggplot2绘制堆积比例图 ----------------------------------------------------------
#数据整理
identical(rownames(result1),group$Sample)
data <- cbind(rownames(result1),result1)
colnames(data)[1] <- "Samples"
data <- melt(data,id.vars = c("Samples"))
colnames(data) <- c('Samples','celltype','proportion')
#开始绘图
mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')
pdf(file="stacked bar chart.pdf",width = 10,height = 8)
ggplot(data,
       aes(Samples,proportion,fill=celltype))+geom_bar(stat="identity",position="fill")+#x 轴是变量 Samples，y 轴是变量 proportion，条形的填充颜色由变量 celltype 决定
       scale_fill_manual(values=mycolors)+#填入需要填充的颜色，22种免疫细胞
       ggtitle("Proportion of immune cells")+theme_gray()+theme(axis.ticks.length=unit(3,'mm'),axis.title.x=element_text(size=11))+
       theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
       guides(fill=guide_legend(title="Types of immune cells"))
dev.off()  


# 盛夏的果实函数科普二 --------------------------------------------------------------
#pivot_longer是tidyr包中的一个函数，它可以将一个数据框中的多个列转换为两列：
#一列是原来的列名，另一列是对应的值。
#pivot_longer(data, cols, names_to = "新的列名", values_to = "新的值列名")

# ggplot2绘制分组箱式图 ----------------------------------------------------------
#数据整理
gp=group
gp=as.data.frame(gp[,2])
row.names(gp)=group[,1]
names(gp)='group'
data1 <- cbind(result1,gp)
colnames(data1)
#data1 <- data1[,c(23,24,1:22)]
data1 <- data1[,c(22,1:21)]
rownames(data1) <- NULL
data1 <- pivot_longer(data = data1,
                      cols = 2:22,
                      names_to = "celltype",
                      values_to = "proportion")
#开始绘图
pdf(file="分组箱式图.pdf",width = 10,height = 8)
ggboxplot(data = data1,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = NULL,#颜色调色板。
          title = "TME Cell composition",#图形标题。
          xlab = NULL,#x 轴标签。
          ylab = "Cell composition",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = FALSE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
dev.off()  


# 提取其中一组（肿瘤/正常）绘图 ---------------------------------------------------------
data2 <- cbind(result1,group)
data2 <- data2[data2$group=='high',]
data2 <- pivot_longer(data = data2,
                      cols = 1:21,
                      names_to = "celltype",
                      values_to = "proportion")
sum(data2$proportion)
pdf(file="单组箱式图.pdf",width = 10,height = 8)
ggboxplot(data = data2,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",
          color = "black",
          xlab = "Types of immune cells",#x 轴标签。
          ylab = NULL,
          title = "TME Cell composition",
          fill = "celltype",
          legend.position = "bottom",
          ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) 
dev.off() 


# ggplot2绘制分组箱式图升级版+显著性P值 -------------------------------------------------
data3 <- cbind(result1,group)
colnames(data3)
data3 <- data3[,c(22,23,1:21)]
data3 <- pivot_longer(data = data3,
                      cols = 3:23,
                      names_to = "celltype",
                      values_to = "proportion")
#开始绘图
pdf(file="分组箱式图+显著性P值.pdf",width = 10,height = 8)
ggboxplot(data = data3,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = c("#81D5B0","#ED5462"),#颜色调色板。
          title = "TME Cell composition",#图形标题。
          xlab = NULL,#x 轴标签。
          ylab = "Cell composition",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = TRUE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
          stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = T) 
dev.off()  
