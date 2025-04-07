
library("stringr")
library("KEGGREST") 
library("data.table")
##################################
gs<-keggGet('hsa04216')
#查找所有基因 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
genelist <- genes[1:length(genes)%%3 ==2] 
genelist <- data.frame(genelist)  

a<-fread("g:/合作文章/WP_FERROPTOSIS.v2023.2.Hs.tsv")
a<-a[17,2]
a<-str_split_fixed(a[1,1],",",65)
a<-a[-1]

b<-c(a,genelist$genelist)
c<-unique(b)

write.csv(c,file = "g:/合作文章/data/Ferroptosis_gene.csv",row.names = F)

a<-fread("g:/合作文章/data/copper_gene.txt",header = F)
b<-str_split_fixed(a[1,1],"、",20)
b<-c(b[1:12],"ATP7B")
write.csv(b,file = "g:/合作文章/data/copper_gene.csv",row.names = F)

a<-fread("g:/合作文章/data/icd_gene.txt",header = F)
b<-str_split_fixed(a[1,1],"、",33)
b<-c(b,"TNF")
write.csv(b,file = "g:/合作文章/data/icd_gene.csv",row.names = F)

a<-read.csv("g:/合作文章/data/Ferroptosis_gene.csv")
a$group<-"Ferroptosis"

b<-read.csv("g:/合作文章/data/copper_gene.csv")
b$group<-"Cuproptosis"

c<-read.csv("g:/合作文章/data/icd_gene.csv")
c$group<-"icd"

d<-rbind(a,b,c)

write.csv(d,file = "g:/合作文章/plot/supplement_table/Supplement Table 4.csv")
