

setwd("/home/data/sdc/wuchx/BioXCG/BioXCG/income/Service/2024/2024.04.28/")
remove(list = ls())
tpm2 <- readRDS("./use_data/BLCA_tpm.rds")


use_tpm <- as.data.frame(t(tpm2[c("LIPT1","ACSL5","CHMP6"),]))
use_tpm$sample <- rownames(use_tpm)
sur_data <- fread("./use_data/survival_data.csv",data.table = F)
use_cox <- dplyr::inner_join(sur_data,use_tpm) %>% na.omit() %>% dplyr::filter(type == "Tumor")
saveRDS(use_cox,"./new_project/use_data/cox_data.rds")


############################## 风险模型构建 ##############################

use_cox <- readRDS("g:/合作文章/data/2024.04.28/new_project/use_data/cox_data.rds")

coxs_model <- use_cox %>% dplyr::select(3:7)
mul_cox_data <- coxs_model
mul_cox <- coxph(Surv(OS.time, OS) ~ ., data = mul_cox_data)
coefficients <- coef(mul_cox)


saveRDS(mul_cox, file = "./new_project/result/model.rds")


# 计算风险评分
model <- readRDS("./new_project/result/model.rds")
model_data <- use_cox %>% dplyr::select(3:7)
model_data$riskScore <- predict(mul_cox, newdata = model_data, type = "risk")

value <- surv_cutpoint(model_data, time = "OS.time", event = "OS", variables = "riskScore") 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])


model_data$Group <- ifelse(model_data$riskScore > cut_off, "High Risk", "Low Risk")
table(model_data$Group)
rownames(model_data) <- use_cox$sample
saveRDS(model_data, file = "./new_project/result/data_exp_sur_riskScore_group.rds")

data_exp_sur_riskScore_group <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")

# 生存分析
data_exp_sur_riskScore_group$Sample <- rownames(data_exp_sur_riskScore_group)
survival_data <- data_exp_sur_riskScore_group[ , c("Sample", "Group", "OS", "OS.time")]
survival_data$OS.time <- round(survival_data$OS.time/30,3)
survival_data$Group <- as.factor(survival_data$Group)

fit <- survfit(Surv(survival_data$OS.time, survival_data$OS) ~ Group, data = survival_data)
fit


# 绘制生存曲线
p <- ggsurvplot(fit,
           pval = TRUE,
           conf.int = F,
           fun = "pct",
           xlab = "Time (Months)",
           palette = c("#bc5148", "#3090a1"),
           legend.title = ggplot2::element_blank(),
           legend.labs = c("High risk Group", "Low risk Group"),
           break.time.by = 20,
           risk.table = T,
           tables.height = 0.2,
           ggtheme = theme_bw())

p
# 添加并居中标题
p$plot <- p$plot + 
  ggtitle("riskScore Group (TCGA-BLCA)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
p



library(timeROC)
data_exp_sur_riskScore_group$OS.time <- data_exp_sur_riskScore_group$OS.time/365
tROC <- timeROC(T = data_exp_sur_riskScore_group$OS.time,
                delta = data_exp_sur_riskScore_group$OS,
                marker = data_exp_sur_riskScore_group$riskScore,
                cause = 1, times = c(1,3,5), ROC = T)
# pdf(file="./result/ROC.pdf", width = 8,height = 6)

plot(tROC, time = 1, col = "red3", title = F, lwd = 2)
plot(tROC, time = 3, col = "green4", title = F, lwd = 2, add = T)
plot(tROC, time = 5, col = "blue3", title = F, lwd = 2, add = T)

legend(0, 1,
       c(paste0("AUC at 1 years  ", round(tROC$AUC[1], 2)),
         paste0("AUC at 3 years  ", round(tROC$AUC[2], 2)),
         paste0("AUC at 5 years  ", round(tROC$AUC[3], 2))),
       col = c("red3", "green4", "blue3"), lwd = 2, cex = 1.2, xjust = -0.9, yjust = 2, bty = "n")



dev.off()

## 基因生存曲线 
data_exp_sur_riskScore_group <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")


data1 <- data_exp_sur_riskScore_group

gene <- "CHMP6"

value <- surv_cutpoint(data1, time = "OS.time", event = "OS", variables = gene) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

single_survival <- data1 %>% 
  dplyr::select(OS,OS.time, gene) %>% ## 选择这四列数据
  dplyr::mutate(group = ifelse(data1[, gene] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>% ## 将时间以月份形式展示
  na.omit() %>% arrange(group) ##按照分组进行排序


single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

# 绘制生存曲线
p <- ggsurvplot(sfit,
           pval = TRUE,
           conf.int = F,
           fun = "pct",
           xlab = "Time (Months)",
           palette = c("#bc5148", "#3090a1"),
           legend.title = ggplot2::element_blank(),
           legend.labs = c(paste0(gene," High"), paste0(gene," Low")),
           break.time.by = 20,
           risk.table = T,
           tables.height = 0.2,
           ggtheme = theme_bw())

p
# 添加并居中标题
p$plot <- p$plot + 
  ggtitle(paste0(gene, " Expression (TCGA-BLCA)")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
p



## 基因表达
exp_tpm <- readRDS("./use_data/BLCA_tpm.rds")


gene <- "ACSL5"
use_data <- as.data.frame(t(exp_tpm[gene,]))
use_data$sample <- rownames(use_data)
use_data <- use_data %>% dplyr::mutate(type = ifelse(grepl("11A$",use_data$sample),"Normal","Tumor"))
table(use_data$type)
library(ggprism)
ggplot(data=use_data,aes(x=type, y=ACSL5, fill=type))+
  geom_boxplot()+
  stat_compare_means(size = 5, comparisons = list(c("Tumor","Normal")),label = "p.format",
                     method = "t.test") +
  ggprism::theme_prism(border = T)+
  labs(y=paste0(gene," Expression (log2(TPM + 1))"),x= NULL,title = "")



use_data <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")

library(ggprism)
gene <- "LIPT1"
ggplot(data=use_data,aes(x=Group, y=LIPT1, fill=Group))+
  geom_boxplot()+
  stat_compare_means(size = 5, comparisons = list(c("Low Risk","High Risk")),label = "p.format",
                     method = "t.test") +
  ggprism::theme_prism(border = T)+
  labs(y=paste0(gene," Expression (log2(TPM + 1))"),x= NULL,title = "")



## 数据集验证

## GSE31684
GSE31684_exp <- readRDS("./result2/GEO/GSE31684_exp.rds")
GSE31684_cli <- readRDS("./result2/GEO/GSE31684_cli.rds")

FADS2 <- as.data.frame(t(GSE31684_exp[c("LIPT1","ACSL5","CHMP6"),]))
FADS2$geo_accession <- rownames(FADS2)

use_data <- dplyr::inner_join(GSE31684_cli,FADS2)
use_data1 <- use_data[,4:8]

model <- readRDS("./new_project/result/model.rds")

model_data <- use_data1
model_data$riskScore <- predict(model, newdata = model_data, type = "risk")


value <- surv_cutpoint(model_data, time = "OS.time", event = "OS", variables = "riskScore") 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

model_data$Group <- ifelse(model_data$riskScore > cut_off, "High Risk", "Low Risk")
table(model_data$Group)
rownames(model_data) <- use_data$geo_accession
saveRDS(model_data, file = "./new_project/result/GSE31684_riskScore_group.rds")




value <- surv_cutpoint(model_data, time = "OS.time", event = "OS", variables = "riskScore") 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

single_survival <- model_data %>% 
  dplyr::select(OS,OS.time, "riskScore") %>% ## 选择这四列数据
  dplyr::mutate(group = ifelse(model_data[, "riskScore"] > cut_off,"High","Low")) %>%
  # mutate(OS.time=round(OS.time/30,2)) %>% ## 将时间以月份形式展示
  na.omit() %>% arrange(group) ##按照分组进行排序


# fwrite(single_survival,"./new_result/GSE_FADS2_group_infor.csv")

single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

# 绘制生存曲线
p <- ggsurvplot(sfit,
                pval = TRUE,
                conf.int = F,
                fun = "pct",
                xlab = "Time (Months)",
                palette = c("#bc5148", "#3090a1"),
                legend.title = ggplot2::element_blank(),
                legend.labs = c("riskScore High", "riskScore Low"),
                break.time.by = 20,
                risk.table = T,
                tables.height = 0.2,
                ggtheme = theme_bw())

p
# 添加并居中标题
p$plot <- p$plot + 
  ggtitle("riskScore (GSE31684)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
p



### 富集分析




### 单基因富集结果

group_infor <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")
tpm_exp <- readRDS("./use_data/BLCA_tpm.rds")
tpm_exp <- tpm_exp %>% dplyr::select(rownames(group_infor))

## 做相关性
use_data <- as.data.frame(t(tpm_exp))
## 创建一个相关性结果的文件
cor_result <- data.frame(matrix(NA,nrow = ncol(use_data),ncol = 3))
colnames(cor_result) <- c("gene","cor","pvalue")
rownames(cor_result) <- colnames(use_data)


all_gene <- c("LIPT1","ACSL5","CHMP6")
# 选择目标基因
use_gene <- "CHMP6"
for (gene_name in colnames(use_data)) {
  
  cor_value <- NA
  p_value <- NA
  
  # 相关性计算
  cor_value <- cor(as.numeric(use_data[,use_gene]), 
                   as.numeric(use_data[,gene_name]),method = "spearman")
  
  # 相关性P值
  p_value <- cor.test(as.numeric(use_data[,use_gene]), 
                      as.numeric(use_data[,gene_name]),method = "spearman")$p.value
  
  # 将结果保存
  cor_result[gene_name,1] <- gene_name
  cor_result[gene_name,2] <- cor_value
  cor_result[gene_name,3] <- p_value
}


# 计算FDR值
cor_result$FDR <- p.adjust(cor_result$pvalue,method = "fdr") 
cor_result1 <- na.omit(cor_result)
# 筛选显著相关的基因（筛选条件根据自己情况调整）
fwrite(cor_result1,"./new_project/result/CHMP6结果/CHMP6_cor_result.csv")




### up基因kegg
remove(list = ls())
# BiocManager::install("GO.db")
library(GO.db)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)



cor_result <- fread("./new_project/result/CHMP6结果/CHMP6_cor_result.csv")
use_gene <- cor_result %>% dplyr::filter(abs(cor) > 0.25 & FDR < 0.05)
#Gene名转化为GeneID
gene.df <- bitr(use_gene$gene, fromType = "SYMBOL",
                toType = c("ENTREZID", "SYMBOL"), 
                OrgDb = org.Hs.eg.db) 


## GO富集分析
GO_all <- clusterProfiler::enrichGO(gene = gene.df$ENTREZID,  
                                    keyType = "ENTREZID",  
                                    OrgDb= org.Hs.eg.db,  
                                    ont = "ALL",   
                                    pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr",  
                                    minGSSize = 10,  
                                    maxGSSize = 500,  
                                    qvalueCutoff = 0.05,  
                                    readable = TRUE)  
GO_result <- data.frame(GO_all)     
table(GO_result$ONTOLOGY)
GO_result <- GO_result %>% dplyr::filter(pvalue < 0.05)
GO_result$GeneRatio <- as.numeric(lapply(GO_result$GeneRatio, function(x) {eval(parse(text = x))}))

# GO三种类别，每种选择显著性最高的12个展示出来
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -pvalue)

fwrite(GO_result,"./new_project/result/CHMP6结果/CHMP6_go_result.csv")



# 纵向柱状图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="CHMP6 cor-gene GO Pathway Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()




## KEGG
kegg <- enrichKEGG(gene.df$ENTREZID, organism = 'human', pvalueCutoff = 0.05, 
                   pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 0.05, 
                   use_internal_data = FALSE)
kegg_result <- data.frame(kegg)

significant_pathways <- kegg_result %>% dplyr::filter(pvalue < 0.05)
fwrite(significant_pathways ,"./new_project/result/CHMP6结果/CHMP6_kegg_result.csv")

# 设置颜色的阶梯和对应的颜色
significant_pathways$pvalue_group <- cut(significant_pathways$p.adjust, breaks = c(0, 0.005, 0.01, 0.05), labels = c("p < 0.005", "0.005 <= p < 0.01", "0.01 <= p < 0.05"))
significant_pathways <- significant_pathways[1:20,]
# 绘制基于ggplot2的气泡图
ggplot(significant_pathways, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue_group), alpha = 1) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "CHMP6 cor-gene KEGG Pathway Enrichment Analysis",x = "Pathway",y = "-log10(P-value)",size = "Count",color = "P-value") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_color_manual(values = c("p < 0.005" = "#DC0000B2", "0.005 <= p < 0.01" = "#F39B7FB2", "0.01 <= p < 0.05" = "#4DBBD5B2"))




## 分组差异

sample_infor <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")
sample_infor$Group <- factor(sample_infor$Group,levels = c("Low Risk", "High Risk"))
sample_infor <- sample_infor %>% dplyr::arrange(Group)
exprSet_by_group <- readRDS("./use_data/BLCA_count.rds")
exprSet_by_group <- exprSet_by_group %>% dplyr::select(rownames(sample_infor))

#edgeR差异分析
group_list <- sample_infor$Group
group_list = factor(group_list)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(exprSet_by_group)
colnames(design) <- levels(group_list)

##差异分析
DGElist <- DGEList(counts = exprSet_by_group, group = group_list)
## 使用cpm值对低表达量的基因进行过滤
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 10 ## 前面做过过滤，这里可做，也可以不做
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
nrDEG_edgeR$gene_id <- rownames(nrDEG_edgeR)
nrDEG_edgeR <- nrDEG_edgeR %>% dplyr::select(gene_id,everything())
#筛选显著性差异的基因
#这里使用logFC > 2 且FDR < 0.05的作为差异基因，可以根据需求改变阈值大小
nrDEG_edgeR_signif <- nrDEG_edgeR %>% filter(abs(logFC) > 1) %>% filter(FDR < 0.05)

fwrite(nrDEG_edgeR,"./new_project/result/风险分组结果/Group_nrDEG_edgeR.csv")
fwrite(nrDEG_edgeR_signif,"./new_project/result/风险分组结果/Group_nrDEG_edgeR_signif.csv")


## 火山图


## LIHC
library(data.table)
nr_DEG <- fread("./new_project/result/风险分组结果/Group_nrDEG_edgeR.csv",data.table = F)


##给差异基因打标签，logFC > 2且 FDR < 0.05认为是上调基因，logFC < -2且 FDR < 0.05认为是下调基因，其它为非差异基因
##这个筛选差异化基因的条件根据自己的情况来定，如果数据量较少可以放宽阈值。
nr_DEG$log10FDR <- -log10(nr_DEG$FDR)
nr_DEG <- nr_DEG %>% 
  mutate(DEG = case_when(logFC > 1 & FDR < 0.05 ~ "up (1006)",
                         abs(logFC) < 1 | FDR > 0.05 ~ "no (20203)",
                         logFC < -1 & FDR < 0.05 ~ "down (493)"))

##加载绘图需要的R包
library(ggplot2) 
library(ggprism) 
library(ggrepel) 
ggplot(nr_DEG, aes(x =logFC, y=log10FDR, colour=DEG)) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c('steelblue','gray','brown')) +
  xlim(-8,8) +  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) + 
  labs(x="log2FC", y="-log10FDR") +
  ggtitle("BLCA DEG") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)



nr_DEG <- fread("./new_project/result/风险分组结果/Group_nrDEG_edgeR_signif.csv",data.table = F)

#Gene名转化为GeneID
gene.df <- bitr(nr_DEG$gene_id, fromType = "SYMBOL",
                toType = c("ENTREZID", "SYMBOL"), 
                OrgDb = org.Hs.eg.db) 


## GO富集分析
GO_all <- clusterProfiler::enrichGO(gene = gene.df$ENTREZID,  
                                    keyType = "ENTREZID",  
                                    OrgDb= org.Hs.eg.db,  
                                    ont = "ALL",   
                                    pvalueCutoff = 0.05, 
                                    pAdjustMethod = "fdr",  
                                    minGSSize = 10,  
                                    maxGSSize = 500,  
                                    qvalueCutoff = 0.05,  
                                    readable = TRUE)  
GO_result <- data.frame(GO_all)     
table(GO_result$ONTOLOGY)
GO_result <- GO_result %>% dplyr::filter(pvalue < 0.05)
GO_result$GeneRatio <- as.numeric(lapply(GO_result$GeneRatio, function(x) {eval(parse(text = x))}))

# GO三种类别，每种选择显著性最高的12个展示出来
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -pvalue)

fwrite(GO_result,"./new_project/result/风险分组结果/分组_go_result.csv")



# 纵向柱状图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="GO Pathway Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()




## KEGG
kegg <- enrichKEGG(gene.df$ENTREZID, organism = 'human', pvalueCutoff = 0.05, 
                   pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 0.05, 
                   use_internal_data = FALSE)
kegg_result <- data.frame(kegg)

significant_pathways <- kegg_result %>% dplyr::filter(pvalue < 0.05) %>% na.omit() 
fwrite(significant_pathways ,"./new_project/result/风险分组结果/分组_kegg_result.csv")

# 设置颜色的阶梯和对应的颜色
significant_pathways$pvalue_group <- cut(significant_pathways$p.adjust, breaks = c(0, 0.005, 0.01, 0.05), labels = c("p < 0.005", "0.005 <= p < 0.01", "0.01 <= p < 0.05"))
# 绘制基于ggplot2的气泡图
ggplot(significant_pathways, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue_group), alpha = 1) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "KEGG Pathway Enrichment Analysis",x = "Pathway",y = "-log10(P-value)",size = "Count",color = "P-value") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_color_manual(values = c("p < 0.005" = "#DC0000B2", "0.005 <= p < 0.01" = "#F39B7FB2", "0.01 <= p < 0.05" = "#4DBBD5B2"))


#### GSEA





## GSEA富集分析
nr_DEG <- fread("./new_project/result/风险分组结果/Group_nrDEG_edgeR.csv",data.table = F)
geneSet <- read.gmt("/home/data/sdc/wuchx/BioXCG/BioXCG/income/Service/2023/2023.11.28/raw_data/h.all.v2023.1.Hs.symbols.gmt") #下载的基因集


# 使用logFC作为排序依据

geneList <- nr_DEG$logFC 
names(geneList) <-  nr_DEG$gene_id
geneList <- sort(geneList, decreasing = T) 

GSEA_enrichment <- GSEA(geneList, 
                        TERM2GENE = geneSet, 
                        pvalueCutoff = 1, 
                        minGSSize = 5, 
                        maxGSSize = 500, 
                        eps = 0, 
                        pAdjustMethod = "BH")

result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)

fwrite(result,"./new_project/result/风险分组结果/分组_GSEA_hallmarker_result.csv")
saveRDS(GSEA_enrichment,"./new_project/result/风险分组结果/分组_GSEA_enrichment.rds")

## 展示最显著的15个通路

# 根据p-value值排序，选择显著的通路或功能（这里设定阈值为0.05）
significant_pathways <- subset(result, pvalue < 0.05)


dotplot(GSEA_enrichment,showCategory=15,color="p.adjust") 

rownames(significant_pathways) <- 1:nrow(significant_pathways)
use_pathway <- significant_pathways[c(1:3,14,21),]


#同时看多个通路
library(enrichplot)
gseaplot2(GSEA_enrichment,c(use_pathway$Description),
          color=brewer.pal(n = 11, name = "RdYlBu")[c(9,10,3,2,1)],pvalue_table = T)



## 免疫浸润

## 免疫浸润等分析
data_exp_sur_riskScore_group <- readRDS("./new_project/result/data_exp_sur_riskScore_group.rds")
data_exp_sur_riskScore_group$sample <- rownames(data_exp_sur_riskScore_group)


immu_data <- fread("./infiltration_estimation_for_tcga.csv",data.table = F)
immu_data$cell_type <- paste0(immu_data$cell_type,"A")


cibersort <- immu_data %>% dplyr::select(cell_type,grep("CIBERSORT$",colnames(immu_data)))
## 删除_CIBERSORT
new_strings <- lapply(colnames(cibersort)[2:23], function(x) gsub("_CIBERSORT", "", x))
new_strings <- unlist(new_strings)
colnames(cibersort) <- c("sample",new_strings)

use_data <- dplyr::inner_join(data_exp_sur_riskScore_group,cibersort)
use_data1 <- use_data[,c(8,7,9:30)] %>% dplyr::arrange(Group)
colnames(use_data1)[2] <- "group"


saveRDS(use_data1,"./new_project/result/风险分组结果/分组_cibersort.rds")

use_data1 <- readRDS("./new_project/result/风险分组结果/分组_cibersort.rds")
use_data2 <- use_data1 %>% tidyr::gather(key = "Variable", value = "Value", -c(sample,group))

library(ggplot2)
library(ggpubr)

ggplot(use_data2, aes(x = Variable, y = Value, fill = group))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test",hide.ns = F)+
  ggprism::theme_prism(border = T)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1) )+
  labs(y="Composition",x= NULL,title = "CIBERSORT")



library(data.table)
library(ggplot2)
library(ggsci)
library(ggprism)
library(scales)

use_data3 <- use_data2 %>% dplyr::arrange(group)

pal<-brewer.pal(12,'Set3')
pal1<-brewer.pal(10,'Paired')


mypal <- append(pal, pal1)
mypal
show_col(mypal,ncol = 5)

ggplot(data=use_data3,aes(sample,Value,fill=Variable))+ 
  geom_bar(stat="identity",position="stack", linewidth=0.75, size=0.25)+
  scale_fill_manual(values = mypal)+
  ggprism::theme_prism(border = T)+
  labs(y="Relative Percent",x = "sample")+
  theme(axis.title.x = element_blank(),  # 移除 x 轴标题
        axis.text.x = element_blank()) 


cor_data <- use_data1[,-c(1,2)]
library(data.table) # 数据出去
library(pheatmap) #热图
library(ComplexHeatmap) #热图
library(RColorBrewer) #颜色
library(circlize) #颜色



use_cor <- cor(cor_data,method = "pearson") #计算相关性
ComplexHeatmap::Heatmap(use_cor,cluster_rows = F,cluster_columns = F,row_names_side = "left",name = "cor",
                        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", use_cor[i, j]), x, y, gp = gpar(fontsize = 7))}
)

colors <- colorRampPalette(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50) ##蓝到红
values <- seq(-1, 1, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)

ComplexHeatmap::Heatmap(use_cor,cluster_rows = F,cluster_columns = F,row_names_side = "left",name = "cor",col = col_fun,
                        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", use_cor[i, j]), x, y, gp = gpar(fontsize = 7))}
)

pheatmap::pheatmap(use_cor, cluster_cols = F, cluster_rows = F, display_numbers = TRUE)


## estimate结果

pancaner_estimate <- readRDS("./pancancer_estimate_score.rds")
pancaner_estimate$sample <- paste0(pancaner_estimate$sample,"A")

estimate <- dplyr::inner_join(data_exp_sur_riskScore_group,pancaner_estimate)
estimate_result <- estimate[,c(8,7,9:11)]

draw_data <-  estimate_result %>% tidyr::gather(key = "Variable", value = "Value", -c(sample,Group))

ggplot(draw_data, aes(x = Variable, y = Value, fill = Group))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test",hide.ns = F)+
  ggprism::theme_prism(border = T)+
  labs(y="Score",x= NULL,title = "ESTIMATE")

## 药物敏感性

drug_data <- readRDS("g:/合作文章/data/2024.04.28/result2/drug_result.rds")
drug_data1 <- drug_data %>% dplyr::select(-Group)
infor <- readRDS("g:/合作文章/data/2024.04.28/new_project/result/data_exp_sur_riskScore_group.rds")
infor$sample <- rownames(infor)
infor <- infor[,7:8]

box_data <- dplyr::inner_join(drug_data1,infor)
library(tidyverse)
risk_group <- box_data %>% pivot_longer(-c(sample, Group), names_to = "Drug", values_to = "Expression")




box_data_high <- box_data %>% filter(box_data$Group == "High Risk")
box_data_low <- box_data %>% filter(box_data$Group == "Low Risk")


# 进行 Wilcoxon 秩和检验
wilcox_test_results <- lapply(1:251, function(i) {
  wilcox.test(box_data_high[,i], box_data_low[,i])
})


# 创建一个空的数据框来存储结果
result_df <- data.frame(
  Drug = character(),
  wilcoxon_test_p_value = numeric(),
  stringsAsFactors = FALSE
)

# 对每个特征进行 t 检验和 Wilcoxon 秩和检验，并将结果存储在数据框中
for (i in 1:length(wilcox_test_results)) {
  result_df <- rbind(result_df, 
                     data.frame(
                       Drug = colnames(box_data)[i],
                       wilcoxon_test_p_value = wilcox_test_results[[i]]$p.value
                     ))
}


p_0.05 <- result_df %>% filter(result_df$wilcoxon_test_p_value < 0.05)
p_0.05$FDR <- p.adjust(p_0.05$wilcoxon_test_p_value)
p_0.05 <- p_0.05 %>% dplyr::filter(FDR < 0.05) %>% dplyr::arrange(FDR)


a <- box_data %>% dplyr::select(Group,everything()) %>% dplyr::select(-sample)
b <- a %>% group_by(Group) %>% summarise_all(mean)
c <- b[,-1]
rownames(c) <- b$Group
d <- as.data.frame(t(c))
d$result <- d$`High Risk` - d$`Low Risk`
d$Drug <- rownames(d)

use_drug <- dplyr::inner_join(p_0.05,d)

# 仅展示显著的药物
# box_data_p <- box_data %>% dplyr::select(p_0.05$Drug)
box_data_p <- box_data %>% dplyr::select(Pyrimethamine,FH535,Methotrexate,Masitinib,Sunitinib,Dasatinib)
box_data_p$sample <- box_data$sample
box_data_p$Group <- box_data$Group

risk_group_p <- box_data_p %>% pivot_longer(-c(sample, Group), names_to = "Drug", values_to = "Expression")
risk_group_p$Drug <- factor(risk_group_p$Drug,levels = c("Pyrimethamine","FH535","Methotrexate","Masitinib","Sunitinib","Dasatinib"))


ggplot(risk_group_p, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() + facet_wrap(~Drug, ncol = 3, scales = "free_y") +
  ylab("Predicted Drug Sensitivity") + xlab("") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "Group") +
  geom_signif(comparisons = list(c("High Risk","Low Risk")),
              map_signif_level = TRUE, test = wilcox.test) +
  cowplot::theme_cowplot(font_size = 15, line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))







## 相关性图

infor <- readRDS("g:/合作文章/data/2024.04.28/new_project/result/data_exp_sur_riskScore_group.rds")
infor$sample <- rownames(infor)
infor <- infor[,c(8,6)]

cor_data <- dplyr::inner_join(infor,box_data_p)
cor_data1<- cor_data %>% dplyr::select(riskScore,Pyrimethamine,FH535,Methotrexate,Masitinib,Sunitinib,Dasatinib)

cor_data2 <- cor_data1 %>% pivot_longer(-c(riskScore), names_to = "Drug")


ggplot(cor_data2, aes(x = value, y = riskScore))+
  geom_point()+
  geom_smooth(method = "lm", fullrange = TRUE) +
  facet_wrap(~Drug)+
  theme_bw()+
  stat_cor(method = "pearson", label.x = 3, label.y = 10,size = 6)+
  theme(strip.text = element_text(size = 15))+
  labs(title = "riskScore Drug Sensitivity",x = "Drug",y = "riskScore")+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 18),    # 调整 x 轴标题字体大小
        axis.title.y = element_text(size = 18),    # 调整 y 轴标题字体大小
        axis.text.x = element_text(size = 14),     # 调整 x 轴刻度标签字体大小
        axis.text.y = element_text(size = 14)      # 调整 y 轴刻度标签字体大小
  )


## 突变分析



# 加载包
library(tidyverse)
# BiocManager::install("maftools")
library(maftools)


# 加载数据
tcga_maf <- data.table::fread("g:/合作文章/data/cohortMAF.2024-06-20.gz")
# `mc3.v0.2.8.PUBLIC.maf` can be downloaded from [https://gdc.cancer.gov/](https://gdc.cancer.gov/).

data_exp_sur_riskScore_group <- readRDS("g:/合作文章/data/2024.04.28/new_project/result/data_exp_sur_riskScore_group.rds")
data_exp_sur_riskScore_group$sample <- rownames(data_exp_sur_riskScore_group)

highrisk <- data_exp_sur_riskScore_group %>% dplyr::filter(Group == "High Risk")
lowrisk <- data_exp_sur_riskScore_group %>% dplyr::filter(Group == "Low Risk")

highrisk_maf <- tcga_maf %>% dplyr::filter(substr(tcga_maf$Tumor_Sample_Barcode, 1, 16) %in% highrisk$sample)
lowrisk_maf <- tcga_maf %>% dplyr::filter(substr(tcga_maf$Tumor_Sample_Barcode, 1, 16) %in% lowrisk$sample)

length(table(highrisk_maf$Tumor_Sample_Barcode))
length(table(lowrisk_maf$Tumor_Sample_Barcode))

highrisk_maf <- read.maf(maf = highrisk_maf)
lowrisk_maf <- read.maf(maf = lowrisk_maf)

# save(highrisk_maf, lowrisk_maf, file = "./data/tmp_data.Rdata")
# load("./data/tmp_data.Rdata")

# 整体展示
plotmafSummary(maf = highrisk_maf, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)
plotmafSummary(maf = lowrisk_maf, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)

# 瀑布图
oncoplot(maf = highrisk_maf)
oncoplot(maf = lowrisk_maf)

# TMB
tmb(maf = highrisk_maf)
tmb(maf = lowrisk_maf)


# laml.titv = titv(maf = highrisk_maf, plot = T, useSyn = T)


# 与TCGA对比
par(mar = c(5, 5, 3, 3)) # 下左上右
tmb_vs_tcga_high = tcgaCompare(maf = highrisk_maf, cohortName = 'High Risk', logscale = TRUE, capture_size = 50)
dev.off()

par(mar = c(5, 5, 3, 3)) # 下左上右
tmb_vs_tcga_low = tcgaCompare(maf = lowrisk_maf, cohortName = 'Low Risk', logscale = TRUE, capture_size = 50)
dev.off()




