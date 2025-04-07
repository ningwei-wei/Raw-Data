library(readxl)
library(tidyverse)

############################ til_score

til_score<-read_xlsx("g:/合作文章/data/1-s2.0-S2211124718304479-mmc2.xlsx")
til_score<-til_score%>%select("ParticipantBarcode","til_percentage")%>%as.data.frame()
colnames(til_score)[1]<-"sample"

TPM<-readRDS("g:/合作文章/data/BLCA_tpm.rds")
TPM<-TPM["Fsp1",]


risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
write.csv(risk_group,file = "g:/合作文章/plot/supplement_table/Supplement Table 2.csv")
risk_group$sample<-substr(rownames(risk_group),1,12)



result<-inner_join(til_score,risk_group,by="sample")

library(ggpubr)
# 设置比较组
my_comparisons <- list( c("High Risk", "Low Risk"))
pdf("g:/合作文章/plot/figure7/til.pdf",width = 8,height = 6)
ggviolin(result, x = "Group", y = "til_percentage", fill = "Group",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
dev.off()



