#绘图
TIDE=read.csv('TIDE.csv',header = T,sep=',')
load(file='E:/R/jiahui_doctor/group.Rdata')
group$sample=row.names(group)
row.names(TIDE)=TIDE[,1]
TIDE=merge(TIDE,group,by='sample')
library(ggplot2)
library(ggpubr)
ggplot(TIDE)+
  geom_histogram(aes(TIDE, fill=group), position = "identity", alpha=0.6, color="white")+
  scale_fill_brewer(palette = "Dark2")

ggplot(TIDE, aes(group, TIDE),color = "group")+
  geom_violin(aes(fill = group))+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "top")+
  stat_compare_means(method = "t.test",label.y = 3)

ggplot(TIDE, aes(TIDE))+
  geom_histogram(aes(y = ..density.., fill = group), position = "identity", alpha = 0.6,color = "white")+
  stat_density(geom = 'line', size=1, position = 'identity', aes(color = group))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  theme_grey(base_size = 22)+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title= element_text(size = 18))+    #修改坐标轴数字大小
  annotate(geom="text", x=-1.3, y=0.75, label="T-test , p = 0.016",size=4.5)+  #图中加入文本
  theme_bw(base_size=12)  #图例大小
#geom_text(x=-2, y=0.6, label="T-test,p = 0.015") 