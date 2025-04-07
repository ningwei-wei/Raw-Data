################### GSVA
library(ggplot2) #绘图使用
library(ComplexHeatmap) #绘图使用
library(clusterProfiler) #数据处理使用
library(GSVA) #GSVA使用
library(GSEABase) #数据处理使用
library(dplyr) #数据处理使用
library(data.table) #数据读取使用

library(survival)
library(survminer)

data_tpm<-readRDS("g:/合作文章/data/BLCA_tpm.rds")
gene_set<-list(c("ABI1", "ARF6", "F11R", "BAIAP2L1", "MYADM", 
            "ACTR3", "ARHGEF5", "SRC", "CAPZA2", "GOLPH3", 
            "RAB10", "ARF4", "CAP1", "PSEN1", "YWHAG", "YWHAZ", "HSPA4"))
gsvapar <- gsvaParam(as.matrix(data_tpm), gene_set, maxDiff=TRUE,kcdf = "Gaussian")

gsva_result <- gsva(gsvapar)

gsva_result<-t(gsva_result)
colnames(gsva_result)<-"TMT"
gsva_result<-as.data.frame(gsva_result)
gsva_result$sample<-rownames(gsva_result)

#risk_group<-read.csv("g:/合作文章/result/riskscore_group.csv")
#colnames(risk_group)[1]<-"sample"

risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
risk_group$sample<-rownames(risk_group)



result<-inner_join(gsva_result,risk_group,by="sample")
############################ 分类
result_height<-result%>%filter(Group=="High Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "TMT",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "TMT") %>%
  dplyr::mutate(group = ifelse(result_height[, "TMT"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
              pval = TRUE,
              conf.int = FALSE,
              fun = "pct",
              xlab = "Time (Months)",
              palette = c("#bc5148", "#3090a1"),
              legend.title = ggplot2::element_blank(),
              legend.labs = c("TMT high", "TMT low"),
              break.time.by = 30,
              risk.table = TRUE,
              tables.height = 0.3,
              ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("High Risk Group") +theme(plot.title = element_text(hjust = 0.5))

pdf("g:/合作文章/plot/figure7/riskscore_height.pdf",width = 8,height = 6)
print(p)
dev.off()
############################################
result_height<-result%>%filter(Group=="Low Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "TMT",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "TMT") %>%
  dplyr::mutate(group = ifelse(result_height[, "TMT"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
              pval = TRUE,
              conf.int = FALSE,
              fun = "pct",
              xlab = "Time (Months)",
              palette = c("#bc5148", "#3090a1"),
              legend.title = ggplot2::element_blank(),
              legend.labs = c("TMT high", "TMT low"),
              break.time.by = 30,
              risk.table = TRUE,
              tables.height = 0.3,
              ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("Low Risk Group") +theme(plot.title = element_text(hjust = 0.5))

pdf("g:/合作文章/plot/figure7/riskscore_low.pdf",width = 8,height = 6)
print(p)
dev.off()

############################### CXCR5+CD8+

data_tpm<-readRDS("g:/合作文章/data/BLCA_tpm.rds")
gene_set<-c("CD8A", "CXCR5", "CXCR3", "ICOS", "CD27", "IL21", "TNF", "TNFRSF6B", "PDCD1",
"TBX21", "SLAMF6","IL27RA")

data_tpm<-data_tpm[gene_set,]
a<-apply(data_tpm, 2, mean)%>%as.data.frame() 
colnames(a)<-"Cxcr5cd8"
a$sample<-rownames(a)

risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
risk_group$sample<-rownames(risk_group)

result<-inner_join(a,risk_group,by="sample")

#########################
result_height<-result%>%filter(Group=="High Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "Cxcr5cd8",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "Cxcr5cd8") %>%
  dplyr::mutate(group = ifelse(result_height[, "Cxcr5cd8"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
           pval = TRUE,
           conf.int = FALSE,
           fun = "pct",
           xlab = "Time (Months)",
           palette = c("#bc5148", "#3090a1"),
           legend.title = ggplot2::element_blank(),
           legend.labs = c("CXCR5+CD8+ high", "CXCR5+CD8+ low"),
           break.time.by = 30,
           risk.table = TRUE,
           tables.height = 0.3,
           ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("High Risk Group") +theme(plot.title = element_text(hjust = 0.5))
pdf(paste0("g:/合作文章/plot/figure7/high.pdf"),width = 8,height = 6)
print(p)
dev.off()
#########################
result_height<-result%>%filter(Group=="Low Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "Cxcr5cd8",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "Cxcr5cd8") %>%
  dplyr::mutate(group = ifelse(result_height[, "Cxcr5cd8"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
              pval = TRUE,
              conf.int = FALSE,
              fun = "pct",
              xlab = "Time (Months)",
              palette = c("#bc5148", "#3090a1"),
              legend.title = ggplot2::element_blank(),
              legend.labs = c("CXCR5+CD8+ high", "CXCR5+CD8+ low"),
              break.time.by = 30,
              risk.table = TRUE,
              tables.height = 0.3,
              ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("Low Risk Group") +theme(plot.title = element_text(hjust = 0.5))
pdf(paste0("g:/合作文章/plot/figure7/low.pdf"),width = 8,height = 6)
print(p)
dev.off()

####################################
library(dplyr)
library(stringr)

TIDE<-read.table("g:/合作文章/data/Results/Tumor_Dysf_Excl_scores/TCGA.BLCA.RNASeq.norm_subtract.OS_base")
TIDE$TIDE<-apply(TIDE[,4:5], 1,max)
TIDE$sample<-substr(rownames(TIDE),1,12)

risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
risk_group$sample<-substr(rownames(risk_group),1,12)


result<-inner_join(TIDE,risk_group,by="sample")

library(ggpubr)

pic <- ggplot(data = result, aes(x = Group, y = TIDE, color = Group))+   #指定数据集，设置坐标轴名称、类别颜色
  geom_boxplot()  +theme_classic() +stat_compare_means()#绘制箱线图

pdf("g:/合作文章/plot/TIDE.pdf",width = 8,height = 6)
pic   #输出箱线图 
dev.off()
#########################
result_height<-result%>%filter(Group=="Low Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "TIDE",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "TIDE") %>%
  dplyr::mutate(group = ifelse(result_height[, "TIDE"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
              pval = TRUE,
              conf.int = FALSE,
              fun = "pct",
              xlab = "Time (Months)",
              palette = c("#bc5148", "#3090a1"),
              legend.title = ggplot2::element_blank(),
              legend.labs = c("TIDE high", "TIDE low"),
              break.time.by = 30,
              risk.table = TRUE,
              tables.height = 0.3,
              ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("Low Risk Group") +theme(plot.title = element_text(hjust = 0.5))

pdf(paste0("g:/合作文章/plot/figure7/TIDE_low.pdf"),width = 8,height = 6)
print(p)
dev.off()

#########################
result_height<-result%>%filter(Group=="High Risk")
value <- surv_cutpoint(result_height, time = "OS.time", event = "OS", variables = "TIDE",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- result_height %>% 
  dplyr::select(OS,OS.time, "TIDE") %>%
  dplyr::mutate(group = ifelse(result_height[, "TIDE"] > cut_off,"High","Low")) %>%
  mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))

#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)

p<-ggsurvplot(sfit,
              pval = TRUE,
              conf.int = FALSE,
              fun = "pct",
              xlab = "Time (Months)",
              palette = c("#bc5148", "#3090a1"),
              legend.title = ggplot2::element_blank(),
              legend.labs = c("TIDE high", "TIDE low"),
              break.time.by = 30,
              risk.table = TRUE,
              tables.height = 0.3,
              ggtheme = theme_bw()) 
p$plot <- p$plot + 
  ggtitle("Hight Risk Group") +theme(plot.title = element_text(hjust = 0.5))

pdf(paste0("g:/合作文章/plot/figure7/TIDE_high.pdf"),width = 8,height = 6)
print(p)
dev.off()



