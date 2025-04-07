library(data.table)
library(tidyverse)
library(stringr)
library(reshape2)
library(ggpubr)
library(ggplot2)

a<-fread("g:/合作文章/data/TCGA-BLCA.GDC_phenotype.tsv.gz",data.table = F)
a<-a[,c(1,5)]
colnames(a)[1]<-"sample"


a$diagnosis_subtype

risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
risk_group$sample<-rownames(risk_group)

b<-left_join(risk_group,a)
b<-table(b$Group,b$additional_treatment_completion_success_outcome)%>%as.data.frame()
for(i in 1:nrow(b)){
  if(b[i,1]=="High Risk"){
    b[i,3]<-b[i,3]/139
  }else{
    b[i,3]<-b[i,3]/261
  }
}
b<-b[-c(1:2),]
p<-ggplot(data=b, mapping=aes(x = Var2, y = Freq,fill=Var1))+
  geom_bar(stat="identity",position=position_dodge(0.9))+theme_classic()+
  xlab("Treatment completion succsecc outcome")+scale_fill_manual(values =c("#E7B800","#4682B4" ))
ggsave(filename = "g:/合作文章/plot/supplement_figure/新建文件夹/Treatment.pdf",width = 8,height = 6)

###################################
a<-fread("g:/合作文章/data/TCGA-BLCA.GDC_phenotype.tsv.gz",data.table = F)
a<-a[,c(1,22)]
colnames(a)[1]<-"sample"


a$diagnosis_subtype

risk_group<-readRDS("g:/合作文章/data/data_exp_sur_riskScore_group.rds")
risk_group$sample<-rownames(risk_group)

b<-left_join(risk_group,a)
b<-table(b$Group,b$diagnosis_subtype)%>%as.data.frame()

for(i in 1:nrow(b)){
  if(b[i,1]=="High Risk"){
    b[i,3]<-b[i,3]/139
  }else{
    b[i,3]<-b[i,3]/261
  }
}
b<-b[-c(1:2),]
p<-ggplot(data=b, mapping=aes(x = Var2, y = Freq,fill=Var1))+
  geom_bar(stat="identity",position=position_dodge(0.9))+theme_classic()+
  xlab("Treatment completion succsecc outcome")+scale_fill_manual(values =c("#E7B800","#4682B4" ))
ggsave(filename = "g:/合作文章/plot/supplement_figure/新建文件夹/Treatment.pdf",width = 8,height = 6)

