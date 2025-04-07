library(UpSetR)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

DEG<-read.csv("g:/合作文章/data/2024.04.28/result/nrDEG_edgeR_signif.csv")
up_gene<-DEG%>%filter(logFC>1,FDR<0.05)%>%select(gene_id)
down_gene<-DEG%>%filter(logFC< -1,FDR<0.05)%>%select(gene_id)


copper_gene<-read.csv("g:/合作文章/data/copper_gene.csv")
ferroto_gene<-read.csv("g:/合作文章/data/Ferroptosis_gene.csv")
icd_gene<-read.csv("g:/合作文章/data/icd_gene.csv")


d1<-list(up_gene=up_gene$gene_id,down_gene=down_gene$gene_id,copper_gene=copper_gene$x)
pdf("g:/合作文章/plot/supplement_figure/copper_gene.pdf",width = 8,height = 6)
upset(fromList(d1),
      number.angles = 0, #交互集合柱状图的柱标倾角
      point.size=4, #图中点的大小
      line.size=1, #图中连接线粗细
      mainbar.y.label="Intersection size", #y轴的标签
      main.bar.color = 'black', #y轴柱状图颜色
      matrix.color="black", #x轴点的颜色
      sets.x.label="Set size",   #x轴的标签
      sets.bar.color=brewer.pal(3,"Set1"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), #bar plot和matrix plot图形高度的占比
      order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale=c(1.5,1.5,1.5,1.5,1.5,2), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color="red" #图中阴影部分的颜色
)
dev.off()

d1<-list(up_gene=up_gene$gene_id,down_gene=down_gene$gene_id,ferroto_gene=ferroto_gene$x)
pdf("g:/合作文章/plot/supplement_figure/ferroto_gene.pdf",width = 8,height = 6)
upset(fromList(d1),
      number.angles = 0, #交互集合柱状图的柱标倾角
      point.size=4, #图中点的大小
      line.size=1, #图中连接线粗细
      mainbar.y.label="Intersection size", #y轴的标签
      main.bar.color = 'black', #y轴柱状图颜色
      matrix.color="black", #x轴点的颜色
      sets.x.label="Set size",   #x轴的标签
      sets.bar.color=brewer.pal(3,"Set1"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), #bar plot和matrix plot图形高度的占比
      order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale=c(1.5,1.5,1.5,1.5,1.5,2), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color="red" #图中阴影部分的颜色
      )
dev.off()

d1<-list(up_gene=up_gene$gene_id,down_gene=down_gene$gene_id,icd_gene=icd_gene$x)
pdf("g:/合作文章/plot/supplement_figure/icd_gene.pdf",width = 8,height = 6)
upset(fromList(d1),
      number.angles = 0, #交互集合柱状图的柱标倾角
      point.size=4, #图中点的大小
      line.size=1, #图中连接线粗细
      mainbar.y.label="Intersection size", #y轴的标签
      main.bar.color = 'black', #y轴柱状图颜色
      matrix.color="black", #x轴点的颜色
      sets.x.label="Set size",   #x轴的标签
      sets.bar.color=brewer.pal(3,"Set1"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      mb.ratio = c(0.7, 0.3), #bar plot和matrix plot图形高度的占比
      order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
      text.scale=c(1.5,1.5,1.5,1.5,1.5,2), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      shade.color="red" #图中阴影部分的颜色
)
dev.off()
      