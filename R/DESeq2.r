###
#* @Date: 2021-08-09 17:26:19
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-08-09 17:27:12
#* @FilePath: /metaSC/R/DESeq2.r
#* @Description:
###
WORK_DIR="/home/pg2020/ug020080910014/Work/2021_Spr-Cources/Omics/lab2-1/"
setwd(WORK_DIR)
#BiocManager::install("DESeq2")
library("DESeq2")
input_data <- read.table("Analyze/counts.txt",header = TRUE,row.names =1,quote="")
input_data <- as.matrix(input_data[,6:9])
condition <- factor(c(rep("time0",2),rep("time1",2)))
coldata <- data.frame(row.names=colnames(input_data),condition)
dds <- DESeqDataSetFromMatrix(countData =input_data,colData=coldata,design= ~condition)
#删除count数小于1或等于1的
dds1 <- dds[ rowSums(counts(dds)) > 1, ]
#进行差异基因分析
dds2 <- DESeq(dds1)
res_0.05 <- results(dds2,alpha=0.05)
res_0.05 <- res_0.05[order(res_0.05$padj),]
resdata_0.05 <- merge(as.data.frame(res_0.05),as.data.frame(counts(dds2,normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata_0.05,file="Analyze/diff_gene_whole.csv",quote=F,row.names=F)
#按照自定义的阈值提取差异基因并导出
diff_gene_deseq2 <-subset(resdata_0.05, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file= "Analyze/DEG.csv",quote=F,row.names=F)

####################################Expression Normalization######################################
#Enter R environment: type "R" and press Enter
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.vM26.annotation.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
Length <- lapply(exons_gene,function(x){sum(width(reduce(x)))}) # time-consuming
Length <-t(as.data.frame(Length))
rownames(Length)<-gsub("\\.(\\.?\\d*)","",rownames(Length))

rownames(input_data)<-gsub("\\.(\\.?\\d*)","",rownames(input_data))
#Genome Version inconsistency
input_data<-as.data.frame(cbind(input_data[match(intersect(rownames(input_data),rownames(Length)),rownames(input_data)),],
                                Length[match(intersect(rownames(input_data),rownames(Length)),rownames(Length)),,drop=F]))
colnames(input_data)[dim(input_data)[2]]<-"Length"

#TPM
kb <- input_data$Length / 1000
countdata <- input_data[,1:(dim(input_data)[2]-1)]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.table(tpm,file="Analyze/Expression_fpkm_tpm.txt",sep="\t",quote=F)
#FPKM
fpkm <- t(t(rpk)/colSums(countdata) * 10^9)
head(fpkm)
write.table(fpkm,file="Analyze/Expression_fpkm.txt",sep="\t",quote=F)
#FPKM to TPM
fpkm_to_tpm = t(t(fpkm)/colSums(fpkm))*10^6
head(fpkm_to_tpm)
