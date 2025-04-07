dir='nw/合作文章/data/' 
samples=list.files( dir ,pattern = 'gz')
samples 
library(data.table)
library(Seurat)
ctList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('.txt.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
  return(ct)
})

#合并样品
lapply(ctList, dim)
tmp =table(unlist(lapply(ctList, rownames)))
cg = names(tmp)[tmp==length(samples)]
bigct = do.call(cbind,
                lapply(ctList,function(ct){ 
                  ct = ct[cg,] 
                  return(ct)
                }))
sce.all=CreateSeuratObject(counts =  bigct, 
                           min.cells = 5,
                           min.features = 300)
sce.all
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident) 
scobj<-sce.all
save(scobj,file = "nw/合作文章/data/scobj.Rdata")
########################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(harmony)
load("nw/合作文章/data/scobj.Rdata")

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^mt-")
png("nw/scRNA/kong/plot/NB/vlnplot_NB2.png")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
scobj <- subset(scobj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA>0)

a<-scobj@meta.data
a<-a[-which(a$orig.ident=="GSM5329919"),]  ###去除正常样本
scobj<-scobj[,rownames(a)]
scobj$orig.ident<-factor(scobj$orig.ident)

scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(scobj),reduction.name = "pca")

scobj <- JackStraw(scobj,dims = 50, num.replicate = 100)
scobj <- ScoreJackStraw(scobj, dims = 1:50)
JackStrawPlot(scobj, dims = 1:50)
ElbowPlot(scobj,ndims = 50,reduction = "pca")

scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:10, reduction.name = "umap")
p1<-DimPlot(scobj, reduction = "umap", label = T)

scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
ElbowPlot(scobj,ndims = 50,reduction = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "umap")
p2<-DimPlot(scobj, reduction = "umap", label = T)

png("nw/合作文章/plot/harmony.png",width = 1000,height = 800)
p1+p2
dev.off()

scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)
scobj <- FindClusters(scobj,resolution = 0.5)

library(clustree)
png("nw/合作文章/plot/cluster.png",width = 1000,height = 800)
clustree(scobj)
dev.off()

scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.5
Idents(scobj) <- "seurat_clusters"

png("nw/合作文章/plot/umap.png",width = 1000,height = 800)
DimPlot(scobj, reduction = "umap", label = T)
dev.off()

png("nw/合作文章/plot/three_gene.png",width = 1200,height = 800)
FeaturePlot(scobj,features = c("ENSG00000176108","ENSG00000197142","ENSG00000144182"),order = T,
            split.by = "orig.ident") ## EPCAM+/KRT8+,上皮细胞
dev.off()
save(scobj,file = "nw/合作文章/data/scobj_umap.Rdata")
load("nw/合作文章/data/scobj_umap.Rdata")

pdf("nw/合作文章/plot/umap.pdf",width = 10,height = 8)
DimPlot(scobj, reduction = "umap", label = T)
dev.off()

FeaturePlot(scobj,features = c("ENSG00000119888","ENSG00000171345")) ## EPCAM+/KRT19+,上皮细胞
FeaturePlot(scobj,features = c("ENSG00000170458","ENSG00000129226")) ## CD14、CD68,骨髓/巨噬细胞
FeaturePlot(scobj,features = c("ENSG00000271503","ENSG00000169442","ENSG00000111796",
                               "ENSG00000167286")) ## CCL5、CD52、KLRB1和CD3D,T细胞
FeaturePlot(scobj,features = c("ENSG00000011465","ENSG00000182578")) ## 成纤维细胞(DCN+/CSF1R+)

scobj.maker<-FindMarkers(scobj,ident.1 = "11",logfc.threshold = 0.25)
FeaturePlot(scobj,features = c("ENSG00000108821","ENSG00000149591")) ## cluster 12;成纤维
FeaturePlot(scobj,features = c("ENSG00000102755","ENSG00000131477")) ## cluster 13;内皮
FeaturePlot(scobj,features = c("ENSG00000167286","ENSG00000169442")) ## cluster 11;T细胞


anno<-c(rep("Epithelial cells",9),"Myeloid cells","Epithelial cells","T cell","Fibroblasts","Endothelial cells")
names(anno)<-levels(scobj)
scobj<-RenameIdents(scobj,anno)

pdf("nw/合作文章/plot/umap_anno.pdf",width = 10,height = 8)
DimPlot(scobj, reduction = "umap", label = T)
dev.off()
pdf("nw/合作文章/plot/three_gene.pdf",width = 15,height = 8)
FeaturePlot(scobj,features = c("ENSG00000144182","ENSG00000197142","ENSG00000176108"),split.by = "orig.ident",
            order = T) ## 
dev.off()
