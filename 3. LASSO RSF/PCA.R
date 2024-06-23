
p=t
pca1 <- prcomp(p[,-ncol(p)],
               center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1) 
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
library(ggplot2)
ggplot(data = df1,aes(x = PC1,y = PC2,color = t$group))+
  #stat_ellipse(aes(fill = t$group),
               #type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 1)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","blue"))+
  scale_colour_manual(values = c("purple","orange"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(p.pca1,filename = "PCA.pdf")



p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=t$risk))+
  stat_ellipse(aes(fill=t$risk),
               type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab1,y=ylab1,color="")+
  guides(fill=F)
p2+scale_fill_manual(values = c("purple","orange"))+
  scale_colour_manual(values = c("purple","orange"))



####pca三维图
library(scatterplot3d)
scatterplot3d(df1[,1:3], # 第1-3主成分
              # 颜色长度要和样本长度一样，且对应！
              #color = rep(c("#00AFBB", "#FC4E07"),each=50),
              color=ifelse(risk_group$group=='high',"#00AFBB","#FC4E07"),
              pch = 20,   ##不同符号
              angle =120,  ##转角度
              lty.hide = 2,
              asp=4  #坐标轴长短比例
)
legend("topleft",c('high','low'),
       fill=c("#00AFBB",  "#FC4E07"),box.col='black')







####相关性玄图
library(tsne)
t=c[,1:8]
#t$group=surv$risk
t$group=res.cat$risk  ###截断值分组

b=t[,1:8]
a<-gsub("\\-","_",colnames(b))
colnames(b)=a
gene_cor <- cor(b, method = 'pearson')
#去除基因的自相关，也就是对角线的值
diag(gene_cor) <- 0
gene_cor  #最终的基因间表达值Pearson相关性矩阵
#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
library(circlize)
chordDiagram(gene_cor, 
             annotationTrack = c('grid', 'name', 'axis'), #绘制外周圆弧区，显示名称和刻度轴
             grid.col = c(CTB_171A8.1 = 'green3', CTD_2371O3.2 = 'red', LINC00240 = 'orange', RP11_126K1.6 = 'purple', 
                          RP11_872J21.3 = 'skyblue', RP3_500L14.2 = 'blue',SNHG10='pink',TLR8_AS1='cyan'
             ), #定义基因颜色
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), #根据相关性大小展示连线的颜色范围
             annotationTrackHeight = c(0.05, 0.1), #名称离圆弧的距离，以及圆弧的宽度
)