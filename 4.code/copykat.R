library(copykat)
library(Seurat)
library(dplyr)


exp.rawdata <- as.matrix(scobj@assays$RNA@counts)

copykat.test <- copykat(rawmat=exp.rawdata, id.type="E", ngene.chr=5, win.size=25, KS.cut=0.05, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=80)
save(copykat.test,file="nw/合作文章/data/copykat.Rdata")

load("nw/合作文章/data/copykat.Rdata")
pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)
#CNA.test<- CNA.test %>% group_by(chrom) %>% do(head(., n = 100))
######################
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

pdf("nw/合作文章/plot/copykat_heatmap.pdf",width = 10,height = 8)
heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =80, method = "euclidean"), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
dev.off()

##################### 对肿瘤细胞再聚类
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.cells <- gsub("-",".",tumor.cells)
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =80, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

png("nw/合作文章/plot/heatmap_tumor.png",width = 1000,height = 800)
heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
dev.off()
###################################
pred.test <- data.frame(copykat.test$prediction)

a<-scobj@meta.data
a$sample<-rownames(a)
colnames(pred.test)[1]<-"sample"

a<-left_join(a,pred.test)
scobj@meta.data$copykat.pred <- a$copykat.pred
scobj@meta.data$copykat.tumor.pred <- rep("normal", nrow(scobj@meta.data))


p1 <- DimPlot(scobj, label = T)
p2 <- DimPlot(scobj, group.by = "copykat.pred")
png("nw/合作文章/plot/copykat.png",width = 1000,height = 800)
p1 + p2 
dev.off()

pdf("nw/合作文章/plot/copykat.pdf",width = 10,height = 8)
DimPlot(scobj, group.by = "copykat.pred")
dev.off()




