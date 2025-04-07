install.packages("BiocManager")
BiocManager::install("hgu133plus2.db")
remove(list = ls()) ##清空当前环境
setwd("C:/Users/Administrator/Desktop/OV") ##设置路径

#1.count、tpm数据整合
library(data.table)
library(dplyr)
library(limma)
library(survival)
##数据处理
##读取样本文件信息
sample_sheet <- fread("./raw_data/gdc_sample_sheet.2024-05-22.tsv")
##取样本ID的前15个字符，例如：TCGA-B6-A0RH-01,作为样本的barcode
sample_sheet$Barcode <- substr(sample_sheet$`Sample ID`,1,15)
##去除重复样本，一个样本的一个组织多次测序的数据
sample_sheet1 <- sample_sheet %>% filter(!duplicated(sample_sheet$Barcode))
##根据Barcode最后两位筛选样本，01表示原发肿瘤样本，11表示正常样本，06表示转移样本
sample_sheet2 <- sample_sheet1 %>% filter(grepl("01$|11$|06$",sample_sheet1$Barcode))
##创建一个数据框，包含"gene_id","gene_name","gene_type"三列，用于后续样本合并
TCGA_OV_Exp <- fread("./raw_data/gdc_download_20240522_080746.843248/00dc8820-694b-47a2-9410-09971711cbee/04eeb9f5-2fd6-4df3-ade8-cdb8d77f4273.rna_seq.augmented_star_gene_counts.tsv")
TCGA_OV_Exp <- TCGA_OV_Exp[!1:4,c("gene_id","gene_name","gene_type")]
##根据样本的文件信息，将所有样本合并成一个数据框COUNTS
for (i in 1:nrow(sample_sheet2)) {
  
  folder_name <- sample_sheet2$`File ID`[i]
  file_name <- sample_sheet2$`File Name`[i]
  sample_name <- sample_sheet2$Barcode[i]
  
  data1 <- fread(paste0("./raw_data/gdc_download_20240522_080746.843248/",folder_name,"/",file_name))
  data2 <- data1[!1:4,c("gene_id","gene_name","gene_type","unstranded")]
  colnames(data2)[4] <- sample_name
  
  TCGA_OV_Exp <- inner_join(TCGA_OV_Exp,data2)
  
}
##数据过滤
##去除50%样本中没有表达的基因，这个过滤条件可以自己筛选，根据个人情况设置
# 计算每行从第四列到最后一列为表达值为0的样本占比
zero_percentage <- rowMeans(TCGA_OV_Exp[, 4:ncol(TCGA_OV_Exp)] == 0)
# 设置阈值为0.5，筛选满足条件基因
TCGA_OV_Exp1 <- TCGA_OV_Exp[zero_percentage < 0.5, ]
#去除第一列和第三列
TCGA_OV_Exp1 <- TCGA_OV_Exp1[ , -c(1, 3)]
#重复的基因取均值
TCGA_OV_Exp2 <- as.data.frame(limma::avereps(TCGA_OV_Exp1[,-1], ID = TCGA_OV_Exp1$gene_name))
##将合并的样本存储起来
saveRDS(TCGA_OV_Exp2,"./use_data/TCGA_OV_counts.rds")

##重新再跑一次，根据样本的文件信息，将所有样本合并成一个数据框TPM
for (i in 1:nrow(sample_sheet2)) {
  
  folder_name <- sample_sheet2$`File ID`[i]
  file_name <- sample_sheet2$`File Name`[i]
  sample_name <- sample_sheet2$Barcode[i]
  
  data1 <- fread(paste0("./raw_data/gdc_download_20240522_080746.843248/",folder_name,"/",file_name))
  data3 <- data1[!1:4,c("gene_id","gene_name","gene_type","tpm_unstranded")]
  colnames(data3)[4] <- sample_name
  
  TCGA_OV_Exp <- inner_join(TCGA_OV_Exp,data3)
  
}
TCGA_OV_tpm <- TCGA_OV_Exp
#TPM标准化log2(x+1)
TCGA_OV_tpm <- as.data.frame(TCGA_OV_tpm)
TCGA_OV_tpm[, 4:ncol(TCGA_OV_tpm)] <- apply(TCGA_OV_tpm[, 4:ncol(TCGA_OV_tpm)], 2, as.numeric)
TCGA_OV_tpm1 <- TCGA_OV_tpm  
TCGA_OV_tpm1[, 4:ncol(TCGA_OV_tpm1)] <- apply(TCGA_OV_tpm1[, 4:ncol(TCGA_OV_tpm1)], 2, function(x) { round(log2(x + 1), 3) })
TCGA_OV_tpm1 <- TCGA_OV_tpm1[, -c(1, 3)]
##数据过滤
##去除50%样本中没有表达的基因，这个过滤条件可以自己筛选，根据个人情况设置
# 计算每行从第二列到最后一列为表达值为0的样本占比
zero_percentage <- rowMeans(TCGA_OV_tpm1[, 2:ncol(TCGA_OV_tpm1)] == 0)
# 设置阈值为0.5，筛选满足条件基因
TCGA_OV_tpm1 <- TCGA_OV_tpm1[zero_percentage < 0.5, ]
#重复的基因取均值
TCGA_OV_tpm2 <- as.data.frame(limma::avereps(TCGA_OV_tpm1[,-1], ID = TCGA_OV_tpm1$gene_name))
saveRDS(TCGA_OV_tpm2,"./use_data/TCGA_OV_tpm.rds")
#-------------------------------------------------------------------------
#临床数据整合
#数据读取
library(dplyr)
library(data.table)
clin_data <- fread("./raw_data/clinical.tsv",data.table = F)
#筛选预后相关的列
use_clin <- clin_data %>% 
  dplyr::select(case_submitter_id,vital_status,days_to_death,days_to_last_follow_up) %>%
  dplyr::filter(!duplicated(case_submitter_id)) ## 去除重复
#把已经去世患者的生存时间和存活患者的生存时间数据合并（生存用0表示，死亡用1表示）
sur_data <- use_clin %>% 
  dplyr::mutate(OS.time = case_when(vital_status == "Alive" ~ days_to_last_follow_up,
                                    vital_status == "Dead" ~ days_to_death)) %>%
  dplyr::mutate(OS = case_when(vital_status == "Alive" ~ 0,
                               vital_status == "Dead" ~ 1))
#去除了含有NA和--的行（去除一些没有生存数据的样本）
sur_data1 <- sur_data %>% select(case_submitter_id,OS,OS.time) %>% na.omit()
sur_data1 <- sur_data1 %>% dplyr::filter(OS.time != "'--")
#保存csv文件
write.csv(sur_data1, file = "./use_data/sur_data1.csv", row.names = TRUE)
#———————————————————————————————————————————————————————————————————————————————
#2.cox分析
#读取
TCGA_OV_tpm1 <- readRDS("./use_data/TCGA_OV_tpm.rds")
sur_data1 <- read.table("./use_data/sur_data1.csv",sep=",",row.names=1,header=T)
use_gene <- fread("./use_data/use_gene.csv")
#把临床信息id加-01
sur_data1$case_submitter_id <- paste0(sur_data1$case_submitter_id, "-01")
#筛选279个基因名所对应的tpm,并于临床信息合并
# 获取数据框的行名
row_names <- rownames(TCGA_OV_tpm1)
# 将行名添加为一个名为 gene_name 的列
draw_data1 <- TCGA_OV_tpm1 %>%
  mutate(gene_name = rownames(TCGA_OV_tpm1)) %>%
  filter(gene_name %in% use_gene$gene_name) %>%
  dplyr::select(-gene_name)
draw_data2 <- as.data.frame(t(draw_data1))
draw_data2$case_submitter_id <- rownames(draw_data2)
draw_data3 <- dplyr::inner_join(sur_data1,draw_data2)
td <- draw_data3[, -1]
#查看类型
#class(td$OS)
#看列名
#colnames(td)
#cox分析
pFilter=1 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("surstat","surtime") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(td[,3:ncol(td)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(OS.time, OS) ~ td[,i], data = td)#开始逐一循环cox分析
  tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  if(pvalue<pFilter){ # 这里我们不需要提取所有基因的数据，只需要有意义的gene的结果，所以设置如果pvalue<0.05才提取数据
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                    cbind(id=i,#合并列，是每个基因的统计数据
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
    )
  }
}
#筛选p<0.1，因为0.05筛选出来太少
outResult1 <- outResult %>% filter(pvalue < 0.1)
#保存
write.table(outResult,file="./use_data/UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
write.csv(outResult1, file = "./use_data/UniCoxSurvival0.1.csv", row.names = F)
#读取
tducs <- read.table("./use_data/UniCoxSurvival0.1.csv",sep=",",row.names=1,header=T)
#提取和格式化p值数据
gene <- rownames(tducs)
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))

pdf(file="./cox_result0.1.pdf", width = 8,height = 6)
#定义图表布局和绘图区域
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))
#绘制基因名称和统计值
xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
#绘制风险比图表
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()
#———————————————————————————————————————————————————————————————————————————————
#3.GSVA富集分析
#数据预处理
library(ggplot2) #绘图使用
library(ComplexHeatmap) #绘图使用
library(clusterProfiler) #数据处理使用
library(GSVA) #GSVA使用
library(GSEABase) #数据处理使用
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
library(BiocParallel)
TCGA_OV_tpm1 <- readRDS("./use_data/TCGA_OV_tpm.rds")
gene_set <- getGmt("./use_data/usegenegmt.gmt")
gsva_data <- TCGA_OV_tpm1
#GSVA富集分析
#创建一个参数对象
params <- gsvaParam(as.matrix(gsva_data), gene_set, minSize = 1, maxSize = Inf, kcdf = "Gaussian", tau = 1, maxDiff = TRUE, absRanking = FALSE)
#执行 GSVA 分析
gsva_result <- gsva(params, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
#转置
gsva_result2 <- as.data.frame(t(gsva_result))
gsva_result2$sample <- rownames(gsva_result2)
#有益基因B-有害基因T
gsva_result2$riskScore <- gsva_result2[,1] - gsva_result2[,2]
write.csv(gsva_result2, file = "./use_data/gsva_result.csv", row.names = F)
#———————————————————————————————————————————————————————————————————————————————
#4.生存分析
#使用最佳cutoff根据riskScore给样本分组
library(survival)
library(survminer)
sur_data1 <- read.table("./use_data/sur_data1.csv",sep=",",row.names=1,header=T)
gsva_result2 <- read.table("./use_data/gsva_result.csv",sep=",",header=T)
#把临床信息id加-01
sur_data1$case_submitter_id <- paste0(sur_data1$case_submitter_id, "-01")
gsva_result2 <- as.data.frame(gsva_result2)
sur_data1 <- as.data.frame(sur_data1)
sur_data2 <- sur_data1 %>%
  dplyr::left_join(gsva_result2 %>% dplyr::select(sample, riskScore), by = c("case_submitter_id" = "sample"))
#选择最佳cutoff
value <- surv_cutpoint(sur_data2, time = "OS.time", event = "OS", variables = "riskScore",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- sur_data2 %>% 
  dplyr::select(case_submitter_id,OS,OS.time, "riskScore") %>%
  dplyr::mutate(group = ifelse(sur_data2[, "riskScore"] > cut_off,"High","Low")) %>%
  #mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)
#pdf("./result/survival_plot.pdf", width = 8, height = 6)
ggsurvplot(sfit,
                pval = TRUE,
                conf.int = FALSE,
                fun = "pct",
                xlab = "Time (Months)",
                palette = "jco",
                legend.title = ggplot2::element_blank(),
                legend.labs = c("Score high", "Score low"),
                break.time.by = 30,
                risk.table = TRUE,
                tables.height = 0.3,
                ggtheme = theme_bw())
#选择group为"low"的行
low_rows <- subset(single_survival, group == "Low")
#选择group为"high"的行
high_rows <- subset(single_survival, group == "High")
#将low的行放在high的行前面
single_survival1 <- rbind(low_rows, high_rows)
write.csv(single_survival1, file = "./use_data/single_survival.csv", row.names = F)
#dev.off()
#-------------------------------------------------------------------------------
#5.差异分析
#数据处理
library(dplyr) # 数据处理
library(data.table) # 数据读取
library(stringr) # 数据处理
library(edgeR) # edgeR差异分析
## 数据读取
TCGA_OV_Exp1 <- readRDS("./use_data/TCGA_OV_counts.rds")
single_survival <- read.table("./use_data/single_survival.csv",sep=",",header=T)
#去除低表达的基因
TCGA_OV_Exp1 = TCGA_OV_Exp1[rowMeans(TCGA_OV_Exp1)>1,] #根据自己的需要去除低表达基因，也可以卡其它阈值
#按照Low和High分组
TCGA_OV_Exp2 <- TCGA_OV_Exp1 %>% dplyr::select(single_survival$case_submitter_id)
#将Low样本和High样本按顺序储存到一个矩阵中
Low_sample <- TCGA_OV_Exp2[,1:204]
High_sample <- TCGA_OV_Exp2[,205:ncol(TCGA_OV_Exp2)]
exprSet_by_group <- cbind(Low_sample,High_sample)
#edgeR差异分析
group_list <- c(rep('Low',ncol(Low_sample)),rep('High',ncol(High_sample)))
##创建分组情况
group_list = factor(group_list)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(exprSet_by_group)
colnames(design) <- levels(group_list)
##差异分析
##DGEList: 这个函数用于创建一个 DGEList 对象，用于存储差异表达分析所需的数据。
##counts 参数用于指定基因的计数数据，group 参数用于指定每个样本所属的组别。
##这样可以将不同组别的基因计数数据整合到一个 DGEList 对象中。
DGElist <- DGEList(counts = exprSet_by_group, group = group_list)
## 使用cpm值对低表达量的基因进行过滤
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 前面做过过滤，这里可做，也可以不做
##将 DGEList 对象中的基因进行子集选择，保留符合条件的基因。
##keep.lib.sizes = FALSE 表示不保留每个样本的库大小信息。
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
##calcNormFactors: 函数用于计算规一化因子（normalization factors）,
##以校正样本间的差异。它根据每个样本中基因的总计数来估计规一化因子。
DGElist <- calcNormFactors( DGElist )
##estimateGLMCommonDisp: 这个函数用于估计共同的离散度（common dispersion），
##它代表所有基因的离散度的平均水平。共同离散度是在样本间共享的离散度估计。
DGElist <- estimateGLMCommonDisp(DGElist, design)
##estimateGLMTrendedDisp: 这个函数用于估计趋势相关的离散度（trended dispersion），
##它表示基因的离散度是否与其平均计数呈现趋势相关。
DGElist <- estimateGLMTrendedDisp(DGElist, design)
##estimateGLMTagwiseDisp: 这个函数用于估计基因特异的离散度（tagwise dispersion），
##它表示每个基因的离散度。基因特异的离散度考虑到了基因间的差异。
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
##这个函数用于拟合广义线性模型，基于之前估计的离散度和设计矩阵。
##它将 DGElist 对象和设计矩阵 design 作为输入，并生成适应于 GLM 的拟合对象 fit
fit <- glmFit(DGElist, design)
##glmLRT: 这个函数用于进行广义线性模型的似然比检验（likelihood ratio test，LRT），
##以获取差异表达的统计显著性。通过指定对比矩阵（contrast），可以比较不同组别之间的差异。
##这里，contrast 设置为 c(-1, 1)，表示比较两个条件之间的差异。
results <- glmLRT(fit, contrast = c(-1, 1))
##topTags: 这个函数用于提取差异表达结果中的前 n 个顶部标记（top tags）。
##在这里，使用 nrow(DGElist) 表示提取全部的顶部标记。结果将存储在 nrDEG_edgeR 中
OV_nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
OV_nrDEG_edgeR <- as.data.frame(OV_nrDEG_edgeR)
write.csv(OV_nrDEG_edgeR, file = "./use_data/OV_nrDEG_edgeR.csv", row.names = T)
#筛选显著性差异的基因
#一般使用logFC > 2 且FDR < 0.05的作为差异基因，这里使用abs(logFC) > 1，不然基因太少
nrDEG_edgeR_signif <- OV_nrDEG_edgeR %>% 
  filter(abs(logFC) > 1) %>% 
  filter(FDR < 0.05)
save(nrDEG_edgeR_signif,file = './use_data/OV_nrDEG_edgeR_signif.Rdata')
#画火山图
library(ggplot2)
library(cowplot)
#读取数据
OV_DEG <- read.table("./use_data/OV_nrDEG_edgeR.csv",sep=",",header=T)
#添加gene_name列
OV_DEG$gene_name <- OV_DEG$X
#删除第一列
OV_DEG <- OV_DEG[, -1]
#把gene_name列移到第一列
OV_DEG <- OV_DEG[, c(ncol(OV_DEG), 1:(ncol(OV_DEG) - 1))]
##给差异基因打标签，logFC > 1且 FDR < 0.05认为是上调基因，logFC < -1且 FDR < 0.05认为是下调基因，其它为非差异基因
##这个筛选差异化基因的条件根据自己的情况来定，如果数据量较少可以放宽阈值。
OV_DEG$log10FDR <- -log10(OV_DEG$FDR)
OV_DEG1 <- OV_DEG %>% 
  mutate(DEG = case_when(logFC > 1 & FDR < 0.05 ~ "up",
                         abs(logFC) < 1 | FDR > 0.05 ~ "no",
                         logFC < -1 & FDR < 0.05 ~ "down"))
write.csv(OV_DEG1, file = "./use_data/OV_DEG.csv", row.names = F)
##随机取500个基因来画
data <- OV_DEG1[sample(nrow(OV_DEG1),500,replace = F),]
data$log10FDR <- data$log10FDR/10 ##调整一下FDR值，纯粹为了图更好看一些
#定义自定义的颜色映射
color_scale <- scale_colour_gradientn(
  colours = c("#0C2C84", "green","yellow", "#CE1256"),
  values = c(0, 0.4,0.7, 1),
  guide = guide_colorbar(title = "-log10FDR"))
#设置点的大小
max_size <- 4
min_size <- 1
data$Size <- sqrt(min_size + (max_size - min_size) * (data$log10FDR - min(data$log10FDR)) / (max(data$log10FDR) - min(data$log10FDR)))
#绘制渐变火山图
ggplot(data , aes(x = logFC, y = log10FDR)) +
  geom_point(aes(color = log10FDR, size = Size)) +
  color_scale +
  labs(x = "Fold Change", y = "-log10(FDR)") +
  ggtitle("OV_DEG") +  # 添加标题
  theme_cowplot() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = c(0.9, 0.8),
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
  geom_text(aes(x = -1, y = 0.7, label = "FC < -1"), vjust = -0.5, hjust = 1.5, size = 3.5) +
  geom_text(aes(x = 1, y = 0.7, label = "FC > 1"), vjust = -0.5, hjust = -0.5, size = 3.5) +
  geom_text(aes(x = 1, y = 0, label = "FDR=0.05"), vjust = 1, hjust = -0.2, size = 3.5) +
  guides(size = "none")+
  xlim(-4, 4)
#-------------------------------------------------------------------------------
#6.富集分析（KEGG\GO\GSEA)
#KEGG富集分析
library(ggplot2) # 绘图使用
library(clusterProfiler) #KEGG富集分析使用
library(org.Hs.eg.db) #转换基因ID使用
library(stats) #数据处理使用
library(data.table) #数据读取使用
library(dplyr) #数据处理使用
#读取数据
OV_DEG <- read.table("./use_data/OV_DEG.csv",sep=",",header=T)
##差异基因筛选,选择abs(logFC) > 2,FDR < 0.05的值作为差异基因。
#这里可以根据自己的情况进行设定abs(logFC) > 1,FDR < 0.05
DEG_data <- OV_DEG %>% dplyr::filter(abs(logFC) > 1 & FDR < 0.05)
DEG_data$gene_id <- DEG_data$gene_name
#Gene名转化为GeneID
gene.df <- bitr(DEG_data$gene_id, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db) #Orgdb是指对应的注释包是哪个
#第一列列名改为gene_id
colnames(gene.df)[1] <- "gene_id"
DEG_data1 <- left_join(gene.df,DEG_data)
#KEGG富集分析（原pvalueCutoff = 0.05和qvalueCutoff = 0.05都改成1，结果太少）
kegg <- enrichKEGG(DEG_data1$ENTREZID, organism = 'human', pvalueCutoff = 1, 
                   pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 1, 
                   use_internal_data = FALSE)

kegg_result <- data.frame(kegg)
kegg_result1 <- kegg_result %>% filter(pvalue < 0.05)
write.csv(kegg_result1, file = "./result/kegg_result.csv", row.names = F)
#画图
kegg_result1 <- read.table("./result/kegg_result.csv",sep=",",header=T)
#条形图
# 根据p-value值排序，选择显著的通路或功能（这里设定阈值为0.05）
significant_pathways <- subset(kegg_result1, pvalue < 0.05)
#选择20个
use_pathway <- significant_pathways[1:20,]
# 绘制基于ggplot2的条形图
ggplot(use_pathway, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Pathway",
       y = "-log10(P-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  cowplot::theme_cowplot()+
  coord_flip()
#气泡图
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(circlize)
library(ggsci)
library(viridis)
use_pathway1 <- use_pathway[,c(4,5,7,11)]
use_pathway1$GeneRatio <- sapply(use_pathway1$GeneRatio, function(x) eval(parse(text = x)))
use_pathway1 <- use_pathway1 %>% dplyr::arrange(GeneRatio)
use_pathway1$Description <- factor(use_pathway1$Description,levels = use_pathway1$Description)
#删去最后一行
use_pathway1 <- use_pathway1[-nrow(use_pathway1), ]

ggplot(use_pathway1, aes(x = GeneRatio, y = Description, size = Count, color = pvalue)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_viridis_c(option = "D") +  # 选择 'D' 作为 Viridis 调色板的选项
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "KEGG Pathway Enrichment Analysis", x = "GeneRatio", y = "Pathway")

#colors <- colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)
#values <- seq(0, 0.05, length.out = 101)[-101]
#col_fun <- colorRamp2(values, colors)
#ggplot(use_pathway1, aes(x = GeneRatio, y = Description, size = Count, color = pvalue)) +
  #geom_point(alpha = 0.9) +
  #scale_size_continuous(range = c(3, 10)) +
  #scale_color_gradientn(colors = col_fun(values)) +
  #theme(plot.title = element_text(hjust = 0.5))+
  #labs(title = "KEGG Pathway Enrichment Analysis", x = "GeneRatio", y = "Pathway")
#-------------------------------------------------------------------------------
#GO富集分析
library(ggpubr) #绘图使用
library(ggplot2) # 绘图使用
library(clusterProfiler) #KEGG富集分析使用
library(org.Hs.eg.db) #转换基因ID使用
library(stats) #数据处理使用
library(data.table) #数据读取使用
library(dplyr) #数据处理使用
#读取数据
OV_DEG <- read.table("./use_data/OV_DEG.csv",sep=",",header=T)
##差异基因筛选,选择abs(logFC) > 2,FDR < 0.05的值作为差异基因。
#这里可以根据自己的情况进行设定abs(logFC) > 1,FDR < 0.05
DEG_data <- OV_DEG %>% dplyr::filter(abs(logFC) > 1 & FDR < 0.05)
DEG_data$gene_id <- DEG_data$gene_name
#Gene名转化为GeneID
gene.df <- bitr(DEG_data$gene_id, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db) #Orgdb是指对应的注释包是哪个
#第一列列名改为gene_id
colnames(gene.df)[1] <- "gene_id"
DEG_data1 <- left_join(gene.df,DEG_data)
#GO富集分析 （pvalueCutoff = 0.01和qvalueCutoff = 0.01改为0.05）
GO_all <- enrichGO(gene = DEG_data1$ENTREZID,  #基因列表(转换的ID)
                   keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                   OrgDb=org.Hs.eg.db,  #物种对应的org包
                   ont = "ALL",   #CC细胞组件，MF分子功能，BF生物学过程，ALL以上三个
                   pvalueCutoff = 0.05,  #p值阈值
                   pAdjustMethod = "fdr",  #多重假设检验校正方式
                   minGSSize = 10,   #注释的最小基因集，默认为10
                   maxGSSize = 500,  #注释的最大基因集，默认为500
                   qvalueCutoff = 0.05,  #q值阈值
                   readable = TRUE)  #基因ID转换为基因名

GO_result <- data.frame(GO_all)  
table(GO_result$ONTOLOGY)
write.csv(GO_result, file = "./result/GO_result.csv", row.names = T)
#画图
#纵向柱状图
#GO三种类别，每种选择显著性最高的12个展示出来
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -p.adjust)
#纵向柱状图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="Top 12 Enriched GO Terms") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()
##横向柱状图
ggplot(go_enrichment_pathway, 
       aes(x=reorder(Description, -Count),y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.6) +  #柱状图宽度
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen")) + #柱状图填充颜色
  facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  labs(x="GO Term", y="Gene_Number", title="Top 12 Enriched GO Terms")+
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 80,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
##绘制气泡图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, Count), y=Count)) +
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_size_continuous(range=c(1, 10)) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  theme_minimal() +
  scale_color_gradient(low = "pink",high ="red")+
  labs(color=expression(-log10(p.adjust),size="Count"), 
       x="Gene Ratio",y="Gene_Number",title="GO Enrichment")+
  theme_bw()
#-------------------------------------------------------------------------------
#GSEA富集分析（根据基因的顺序来计算富集得分的）
library(ggplot2) #画图使用
library(clusterProfiler) #GSEA富集/数据读取使用
library(GSEABase) #GSEA富集使用
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
#读取数据
geneSet <- read.gmt("./raw_data/h.all.v2023.2.Hs.symbols.gmt") #下载的基因集
OV_DEG <- read.table("./use_data/OV_nrDEG_edgeR.csv",sep=",",header=T)
#添加gene_name列
OV_DEG$gene_name <- OV_DEG$X
#删除第一列
OV_DEG <- OV_DEG[, -1]
#把gene_name列移到第一列
OV_DEG <- OV_DEG[, c(ncol(OV_DEG), 1:(ncol(OV_DEG) - 1))]
#直接使用logFC作为排序依据(pvalueCutoff = 0.05改为1)
geneList <- OV_DEG$logFC #获取GeneList
names(geneList) <- OV_DEG$gene_name #使用转换好的gene_name，对GeneList命名
geneList <- sort(geneList, decreasing = T) #从高到低排序
GSEA_enrichment <- GSEA(geneList, #排序后的gene
                        TERM2GENE = geneSet, #基因集
                        pvalueCutoff = 1, #P值阈值（保留所有通路）
                        minGSSize = 10, #最小基因数量
                        maxGSSize = 500, #最大基因数量
                        eps = 0, #P值边界
                        pAdjustMethod = "BH") #校正P值的计算方法

GSEA_result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)
saveRDS(GSEA_enrichment,"./result/GSEA_enrichment.rds")
#根据p-value值排序，选择显著的通路或功能(这里设定阈值为0.05)
significant_pathways <- subset(GSEA_result, pvalue< 0.05)
write.csv(significant_pathways, file = "./result/GSEA_result.csv", row.names = F)
#画图(NES条形图)
library(tidyverse)
library(ggplot2)
library(cols4all)
library(patchwork)
#GSEA_enrichment <- readRDS("./result/GSEA_enrichment.rds")
GSEA_result <- read.table("./result/GSEA_result.csv",sep=",",header=T)
# 添加上下调分组标签
GSEA_result$group <- factor(case_when(
  GSEA_result$NES > 0 ~ 'up',
  GSEA_result$NES < 0 ~ 'down'
), levels = c('up', 'down'))  # 设定分组因子的水平顺序
# 根据分组因子对数据框重新排序，确保"up"组在前面
GSEA_result <- arrange(GSEA_result, group)
# 基础版上下调富集柱形图绘制
p <- ggplot(GSEA_result, aes(x = NES, y =reorder(Description, NES), fill = group)) + 
  geom_col() +  # 绘制条形图
  scale_fill_manual(values = c("up" = "#EE6677", "down" = "#4477AA")) +  # 设置分组颜色
  theme_bw()  # 设置主题
#自定义主题调整：
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60',size = 1.1),
  axis.text = element_text(size = 12)
)
p1 <- p + mytheme
#先根据上下调标签拆分数据框：
up <- GSEA_result[which(GSEA_result$NES > 0),]
down <- GSEA_result[which(GSEA_result$NES < 0),]
#添加下调pathway标签：
p2 <- p1 +
  geom_text(data = down,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 0) #标签右对齐
#添加上调pathway标签：
p3 <- p2 +
  geom_text(data = up,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 1) #标签左对齐
#继续调整细节：
p4 <- p3 +
  scale_x_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 3)) + #x轴刻度修改
  labs(x = 'NES', y = ' ', title = 'Enriched Pathway') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整
#最后添加上下调提示标签：
p5 <- p4 +
  geom_text(x = 2.2, y = 16, label = "Up", size = 6, color = '#EE6677') +
  geom_text(x = -2.5, y = 7, label = "Down", size = 6, color = '#4477AA')
#-------------------------------------------------------------------------------
#7.免疫浸润CIBERSORT分析
library(devtools)
library(preprocessCore)
devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
TCGA_OV_tpm1 <- readRDS("./use_data/TCGA_OV_tpm.rds")
TCGA_OV_tpm1 <- as.matrix(TCGA_OV_tpm1)
data(LM22)
results <- cibersort(sig_matrix = LM22, TCGA_OV_tpm1)
results <- as.matrix(results)
write.csv(results, file = "./result/CIBERSORT_result.csv", row.names = T)
#画图
library(dplyr)
library(ggplot2)
library(ggpubr)
results <- read.table("./result/CIBERSORT_result.csv",sep=",",header=T)
single_survival <- read.table("./use_data/single_survival.csv",sep=",",header=T)
#按照case_submitter_id进行匹配
results2 <- left_join(results, single_survival, by = c("X" = "case_submitter_id"))
results2 <- results2[, -c(24:29)]
#把Low放在前面High放在后面
results2 <- results2[order(results2$group, decreasing = TRUE), ]
#将Low样本和High样本按顺序储存到一个矩阵中
Low_sample <- results2[1:204,]
High_sample <- results2[205:419,]
exprSet_by_group1 <- rbind(Low_sample, High_sample)
use_data2 <- exprSet_by_group1 %>% tidyr::gather(key = "Variable", value = "Value", -c(X,group))
ggplot(use_data2, aes(x = Variable, y = Value, fill = group))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test",hide.ns = F)+
  ggprism::theme_prism(border = T)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1) )+
  labs(y="Composition",x= NULL,title = "CIBERSORT")
#-------------------------------------------------------------------------------
#免疫浸润ESTIMATE分析
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos = rforge, dependencies = TRUE)
# library(estimate)
# library(writexl)
# library(readxl)
# # 将表达矩阵文件转换成为gct文件
# #TCGA_OV_tpm1 <- readRDS("./use_data/TCGA_OV_tpm.rds")
# #write.table(TCGA_OV_tpm1, file = "./use_data/TCGA_OV_tpm.txt", sep = "\t", row.names = T, col.names = TRUE, quote = F)
# exp.file <- read.table("./use_data/TCGA_OV_tpm.txt", header = TRUE, sep = "\t", check.names = FALSE)
# filterCommonGenes(input.f="./use_data/TCGA_OV_tpm.txt", output.f="ov_gene.gct", id="GeneSymbol")
# estimateScore(input.ds = "ov_gene.gct", output.ds = "ov_estimate_score.gct", platform = "illumina")
# #删除前两行。
# scores <- read.table("ov_estimate_score.gct",skip = 2,header = T)
# rownames(scores) <- scores[,1]
# scores <- t(scores[,3:ncol(scores)])
# rownames(scores) <- gsub("\\.", "-", rownames(scores))
# write.csv(scores, file = "./result/ESTIMATE_result.csv", row.names = T)
# #画箱线图
# library(ggplot2)
# library(ggpubr)
# library(dplyr)
# ESTIMATE_result <- read.table("./result/ESTIMATE_result.csv",sep=",",header=T)
# single_survival <- read.table("./use_data/single_survival.csv",sep=",",header=T)
# #按照case_submitter_id进行匹配
# ESTIMATE_result2 <- left_join(ESTIMATE_result, single_survival, by = c("X" = "case_submitter_id"))
# ESTIMATE_result2 <- ESTIMATE_result2[, -c(5:7)]
# #把Low放在前面High放在后面
# ESTIMATE_result2 <- ESTIMATE_result2[order(ESTIMATE_result2$group, decreasing = TRUE), ]
# #将Low样本和High样本按顺序储存到一个矩阵中
# Low_sample1 <- ESTIMATE_result2[1:204,]
# High_sample1 <- ESTIMATE_result2[205:419,]
# exprSet_by_group2 <- rbind(Low_sample1, High_sample1)
# use_data3 <- exprSet_by_group2 %>% tidyr::gather(key = "Variable", value = "Value", -c(X,group))
# #method = "wilcox.test"改成t.test，没有改变，免疫得分效果不好不用此结果
# p7 <- ggplot(use_data3, aes(x = Variable, y = Value, fill = group))+
#   geom_boxplot()+
#   stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test",hide.ns = F)+
#   ggprism::theme_prism(border = T)+
#   theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1) )+
#   labs(y="Composition",x= NULL,title = "ESTIMATE")
#-------------------------------------------------------------------------------
#8.突变数据变异分析
#数据合并
library(R.utils)
library(data.table)
library(dplyr)
##读取样本文件信息
sample_sheet <- fread("./raw_data/MAF/gdc_sample_sheet.2024-05-29.tsv")
tcga_ov_maf <- fread("./raw_data/MAF/gdc_download_20240529_030142.999834/00be7f02-1b20-4dfa-94ee-135d4a092cc1/d0e54c97-b4ec-480e-b1dc-30b8459cf219.wxs.aliquot_ensemble_masked.maf.gz")
tcga_ov_maf <- tcga_ov_maf[1,]
tcga_ov_maf[1,] <- NA

for (i in 1:nrow(sample_sheet)) {
  
  folder_name <- sample_sheet$`File ID`[i]
  file_name <- sample_sheet$`File Name`[i]
  data1 <- fread(paste0("./raw_data/MAF/gdc_download_20240529_030142.999834/",folder_name,"/",file_name),fill = F)
  if (ncol(data1) != 2 ) { # 只有列数不等于2 的数据才进行合并.因为有个别样本没有任何突变，这种样本在读取时，并不是传统140列的格式，而是将头文件读入了。
    tcga_ov_maf <- rbind(tcga_ov_maf,data1)
  }
  
}
tcga_ov_maf <- tcga_ov_maf[-1,]
length(table(tcga_ov_maf$Tumor_Sample_Barcode))
fwrite(tcga_ov_maf,"./use_data/tcga_ov_maf.maf")
#maftools工具分析
library(maftools)
##加载数据
tcga_maf <- data.table::fread("./use_data/tcga_ov_maf.maf")
single_survival <- read.table("./use_data/single_survival.csv",sep=",",header=T)
single_survival <- single_survival[, -c(2:4)]
Low_sample <- single_survival[1:204,]
High_sample <- single_survival[205:419,]
highrisk_maf <- tcga_maf %>% filter(substr(tcga_maf$Tumor_Sample_Barcode, 1, 15) %in% High_sample$case_submitter_id)
lowrisk_maf <- tcga_maf %>% filter(substr(tcga_maf$Tumor_Sample_Barcode, 1, 15) %in% Low_sample$case_submitter_id)
length(table(highrisk_maf$Tumor_Sample_Barcode))
length(table(lowrisk_maf$Tumor_Sample_Barcode))
highrisk_maf <- read.maf(maf = highrisk_maf)
lowrisk_maf <- read.maf(maf = lowrisk_maf)
save(highrisk_maf, lowrisk_maf, file = "./use_data/maf_data.Rdata")
#load("./use_data/maf_data.Rdata")
# 整体展示
plotmafSummary(maf = highrisk_maf, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)
plotmafSummary(maf = lowrisk_maf, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)
# 瀑布图
oncoplot(maf = highrisk_maf,titleText = "OV High Expression")
oncoplot(maf = lowrisk_maf,titleText = "OV Low Expression")
#-------------------------------------------------------------------------------
