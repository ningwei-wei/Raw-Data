
library(ComplexHeatmap)
data_tpm<-readRDS("g:/合作文章/data/BLCA_tpm.rds")
DEG<-read.csv("g:/合作文章/data/2024.04.28/result/nrDEG_edgeR_signif.csv")

all_gene<-read.csv("g:/合作文章/data/all_gene.csv")

copper_gene<-read.csv("g:/合作文章/data/copper_gene.csv")
ferroto_gene<-read.csv("g:/合作文章/data/Ferroptosis_gene.csv")
icd_gene<-read.csv("g:/合作文章/data/icd_gene.csv")

a1<-intersect(ferroto_gene$x,DEG$gene_id)
a2<-intersect(icd_gene$x,DEG$gene_id)


data_tpm<-data_tpm[c(a1,a2),]

heatmap_data<-data_tpm
heatmap_data1 <- as.data.frame(t(heatmap_data))

norm_data <- t(apply(heatmap_data, 1, function(x){(x-mean(x))/sd(x)}))

norm_data[norm_data > 2] <- 2
norm_data[norm_data < -2] <- -2


library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
values <- seq(-2.5, 2.5, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

# 创建注释对象
library(ComplexHeatmap)
top_annotation <- HeatmapAnnotation(Type = c(rep("Normal",19),rep("Tumor",405)), col = list(Type = c("Tumor" = "red","Normal"="blue")))
right_annotation<-HeatmapAnnotation(Type = c(rep("Ferroto_gene",7),rep("Icd_gene",4)), col = list(Type = c("Ferroto_gene" = "red","Icd_gene"="blue")))
# split_vector <- c(rep("", 5), rep(" ",3))
pdf("g:/合作文章/plot/Figure1/heatmap.pdf",width = 10,height = 8)
ComplexHeatmap::Heatmap(norm_data,cluster_rows = F,top_annotation = top_annotation,show_row_dend = F,show_column_dend = F,
                        column_title = "Gene",row_names_side = "left",
                        cluster_columns = F,show_column_names = F,show_row_names = T,row_split = c(rep("Ferroto_gene",7),rep("Icd_gene",4)),
                        show_heatmap_legend = T,name = " ",col = col_fun, column_split = c(rep("Normal",19),rep("Tumor",405)))
dev.off()
