
###############
copper_gene<-read.csv("g:/合作文章/data/copper_gene.csv")
ferroto_gene<-read.csv("g:/合作文章/data/Ferroptosis_gene.csv")
icd_gene<-read.csv("g:/合作文章/data/icd_gene.csv")

all_gene<-rbind(copper_gene,ferroto_gene,icd_gene)
write.csv(a,file="g:/合作文章/data/all_gene.csv",row.names = F,quote = F)

##################################cox分析
library(stringr)
library(dplyr)
library(survival)
data_tpm<-readRDS("g:/合作文章/data/BLCA_tpm.rds")
a<-str_split_fixed(colnames(data_tpm),"-",4)[,4]
b<-which(a=="11A")
data_tpm<-data_tpm[,-b]
data_tpm<-data_tpm[all_gene$x,]
data_tpm<-t(data_tpm)
data_tpm<-as.data.frame(data_tpm)
data_tpm$sample<-rownames(data_tpm)

survi_data<-read.csv("g:/合作文章/data/survival_data.csv")

data<-inner_join(survi_data,data_tpm)
data<-na.omit(data)
rownames(data)<-data$sample
data<-data[,-c(1:2)]


pFilter=0.05 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("surstat","surtime") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(data[,3:ncol(data)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(OS.time, OS) ~ data[,i], data = data)#开始逐一循环cox分析
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
outResult$HR<-as.numeric(outResult$HR)
outResult$L95CI<-as.numeric(outResult$L95CI)
outResult$H95CI<-as.numeric(outResult$H95CI)
outResult$pvalue<-as.numeric(outResult$pvalue)

write.csv(outResult, file = "g:/合作文章/result/single_cox.csv", row.names = F)
outResult<-read.csv("g:/合作文章/result/single_cox.csv")
#提取和格式化p值数据
gene <- outResult$id
hr <- sprintf("%.3f",outResult$"HR")
hrLow  <- sprintf("%.3f",outResult$"L95CI")
hrHigh <- sprintf("%.3f",outResult$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(outResult$pvalue<0.001, "<0.001", sprintf("%.3f", outResult$pvalue))

pdf(file="g:/合作文章/plot/Figure1/single_cox.pdf", width = 10,height = 6)
#定义图表布局和绘图区域
n <- nrow(outResult)
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

################# 多因素cox
formu<-paste0("Surv(OS.time, OS) ~ ",gene[1])
for(i in 2:length(gene)){
  formu<-paste0(formu,"+",gene[i])
  
}
formu<-formula(formu)

outResult<-data.frame()
tdcox <- coxph(formu, data = data)#开始逐一循环cox分析
tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                  cbind(id=rownames(tdcoxSummary$coefficients),#合并列，是每个基因的统计数据
                        HR=tdcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                        L95CI=tdcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                        H95CI=tdcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                        pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
  )
outResult$HR<-as.numeric(outResult$HR)
outResult$L95CI<-as.numeric(outResult$L95CI)
outResult$H95CI<-as.numeric(outResult$H95CI)
outResult$pvalue<-as.numeric(outResult$pvalue)
outResult<-outResult%>%filter(pvalue<0.05)

write.csv(outResult, file = "g:/合作文章/result/multi_cox.csv", row.names = F)
outResult<-read.csv("g:/合作文章/result/multi_cox.csv")
#提取和格式化p值数据
gene <- outResult$id
hr <- sprintf("%.3f",outResult$"HR")
hrLow  <- sprintf("%.3f",outResult$"L95CI")
hrHigh <- sprintf("%.3f",outResult$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(outResult$pvalue<0.001, "<0.001", sprintf("%.3f", outResult$pvalue))

pdf(file="g:/合作文章/plot/Figure1/multi_cox.pdf", width = 10,height = 3)
#定义图表布局和绘图区域
n <- nrow(outResult)
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

############################### KM
library(survival)
library(survminer)

data$riskScore<- -0.64736*data$LIPT1+ -0.19148*data$ACSL5+0.43915*data$CHMP6

#选择最佳cutoff
value <- surv_cutpoint(data, time = "OS.time", event = "OS", variables = "riskScore",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])
#分组
single_survival <- data %>% 
  dplyr::select(OS,OS.time, "riskScore") %>%
  dplyr::mutate(group = ifelse(data[, "riskScore"] > cut_off,"High","Low")) %>%
  #mutate(OS.time=round(OS.time/30,2)) %>%
  na.omit() %>% arrange(group)
single_survival$group <- factor(single_survival$group,levels = c("High","Low"))
single_survival$OS.time<-single_survival$OS.time/30
#绘图
sfit <- survfit(Surv(OS.time, OS) ~ group, data = single_survival)
pdf("g:/合作文章/plot/Figure1/KM.pdf", width = 8, height = 6)
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
dev.off()

write.csv(single_survival, file = "g:/合作文章/result/riskscore_group.csv")

######################### GSE31684验证

gse_data<-fread("g:/合作文章/data/GSE135337_series_matrix.txt.gz",data.table = F,fill = T)

