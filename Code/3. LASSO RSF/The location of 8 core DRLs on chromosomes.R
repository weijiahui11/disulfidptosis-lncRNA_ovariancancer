#install.packages('RCircos')
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram,#这是上面load的基因组文件
                            chr.exclude<- NULL,#这个参数是要排除的染色体，这里我选无，也就是画出所有的染色体，你也可以排除x和Y染色体
                            tracks.inside=10,#这一个参数是指定在染色体圆圈内部一共要画几个圆圈
                            tracks.outside=0)#这个参数指定在外部画几个圆圈
RCircos.Set.Plot.Area()#建立一个画板
RCircos.Chromosome.Ideogram.Plot()

load('E:/R/jiahui_doctor/卵巢癌处理/8lnc染色体位置.RData')
#data=read.csv('CNV染色体位置.csv',header=T,sep=',')
data=lnc8
data=data[,2:5]
data=data[,c(2,3,4,1)]
names(data)[1:4]=c('Chromosome','chromStart','chromEnd','Gene')

name.col <- 4 #数据是4列
side <- "in" #画在基因组骨架的内侧
track.num <- 1#基因组骨架内侧的第一个track位置上画图
RCircos.Gene.Connector.Plot(data,track.num, side)#画connector（连接基因名称和基因组位置）
track.num <- 2
RCircos.Gene.Name.Plot(data, name.col,track.num, side)#加基因名称
