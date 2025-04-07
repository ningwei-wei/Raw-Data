library(data.table)
library(tidyverse)
library(survival)
library(survminer)

All_data<-fread("G:/drug/data/GDC-PANCAN.htseq_fpkm-uq.tsv.gz",data.table = F)
a<-All_data[,-1]
a<-(2^a)-1
All_data[,-1]<-a

FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
a <- apply(a,2,FPKM2TPM)

a1<-grep("ENSG00000144182",All_data$xena_sample)  #### LIPT1
a2<-grep("ENSG00000197142",All_data$xena_sample) ### ACSL5
a3<-grep("ENSG00000176108",All_data$xena_sample) ### CHMP6
All_data<-a[c(a1,a2,a3),]

save(All_data,file = "g:/合作文章/data/pancancer.Rdata")

load("g:/合作文章/data/pancancer.Rdata")

pheno_data<-fread("g:/drug/data/GDC-PANCAN.basic_phenotype.tsv.gz",data.table = F)
pheno_data<-pheno_data%>%filter(sample_type_id==1&project_id!="")

survi_data<-fread("g:/drug/data/GDC-PANCAN.survival.tsv",data.table = F)

pheno_data<-inner_join(pheno_data,survi_data)
pheno_data<-pheno_data%>%filter(sample%in%colnames(All_data))

use_cox <- readRDS("g:/合作文章/data/2024.04.28/new_project/use_data/cox_data.rds")
coxs_model <- use_cox %>% dplyr::select(3:7)
mul_cox_data <- coxs_model
mul_cox <- coxph(Surv(OS.time, OS) ~ ., data = mul_cox_data)


a<-unique(pheno_data$project_id)

for(i in 1:length(a)){
  tmp_pheno<-pheno_data%>%filter(project_id==a[i])
  tmp_data<-All_data[,tmp_pheno$sample]
  tmp_data<-t(tmp_data)%>%as.data.frame()
  colnames(tmp_data)<-c("LIPT1","ACSL5","CHMP6")
  tmp_data$riskScore <-  predict(mul_cox, newdata = tmp_data, type = "risk")
  tmp_data$sample<-rownames(tmp_data)
  tmp_data<-inner_join(tmp_data,tmp_pheno)
  value <- surv_cutpoint(tmp_data, time = "OS.time", event = "OS", variables = "riskScore") 
  cut_off <- as.numeric(value[["cutpoint"]][1, 1])
  single_survival <- tmp_data %>% 
    dplyr::select(OS,OS.time, "riskScore") %>% ## 选择这四列数据
    dplyr::mutate(group = ifelse(tmp_data[, "riskScore"] > cut_off,"High","Low")) %>%
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
                  legend.labs = c("riskScore High", "riskScore Low"),
                  break.time.by = 20,
                  risk.table = T,
                  tables.height = 0.2,
                  ggtheme = theme_bw())
  
  p
  # 添加并居中标题
  p$plot <- p$plot + 
    ggtitle(a[i]) +theme(plot.title = element_text(hjust = 0.5))
  pdf(paste0("g:/合作文章/plot/pancancer/",a[i],".pdf"),width = 8,height = 6)
  print(p)
  dev.off()
}
