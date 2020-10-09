## 0. 安装
# rm(list=ls()) 和rm() 的区别在于后者只能清除指定内容，而前者可以删除所有
rm(list = ls())
options()$repos
options()$BioC_mirror
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos
options()$BioC_mirror

# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("KEGG.db",ask = F,update = F)
BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
# clusterProfiler fault
BiocManager::install("clusterProfiler")
#Error: package or namespace load failed for ‘data.table’ in library.dynam(lib, package, package.lib):
#  shared object ‘datatable.so’ not found
#Error: loading failed
#Execution halted
# install old version.

BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)

# 下面代码被我注释了，意思是这些代码不需要运行，因为它过时了，很多旧教程就忽略
# 在代码前面加上 # 这个符号，代码代码被注释，意思是不会被运行
# source("https://bioconductor.org/biocLite.R")
# library('BiocInstaller')
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# BiocInstaller::biocLite("GEOquery")
# BiocInstaller::biocLite(c("limma"))
# BiocInstaller::biocLite(c("impute"))

# 但是接下来的代码又需要运行啦
options()$repos
install.packages('WGCNA')
# package ‘preprocessCore’ is not available (for R version 3.6.1)
install.packages(c("FactoMineR", "factoextra"))
install.packages(c("ggplot2", "pheatmap","ggpubr"))
library("FactoMineR")
library("factoextra")

library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)


## 1. T 检验
getwd()
df <- read.csv("/home/hwrn/WorkSpace/Qiandaohu_MG/Analyze/ko_tpm.csv", header = T)
head(df)
dim(df)
x1902<-df$X1902
x21<-df$X21
# t test for different in mean value between two samples
# if p-value < 0.05, then aaccept the `alternative hypothesis`
t.test(x1902, x21)
t.test(x1902, x21, alternative = "greater")
t.test(x1902, x21, alternative = "less")
t.test(x1902, x21, alternative = "greater", var.equal = T)


## 2. X^2 test
# zero: rows and cols are not connected
# if p <0.05, select the other assumption
x<-c(60,3,32,11)
x<-c(41, 30, 47, 67)
df <- matrix(x, ncol=2)
chisq.test(df)

## 3. 单因素方差分析
df <- read.csv("./example_data/anova_example.csv", header = T, sep = "\t")
# for using aov(), you should use reshape2 to transform "wide format" to "long format"
library("reshape2")
df1<-melt(df)
one.way<-aov(value~variable, data=df1)
summary(one.way)
# if `Pr(>F)` < 0.05: there are differences, but which? unknown
TukeyHSD(one.way)
# if `p adj` < 0.05, then the two is different


## 4. 双因素分析
df <- read.csv("example_data/two_way_anova.csv", header = T)
head(df)
dim(df)
two.way<-aov(yield~density+fertilizer, data = df)
summary(two.way)
# Then, consider coeffection between two vars like "density" and "fertilizer"
# to do this, change `+` to `*`
two.way1<-aov(yield~density*fertilizer, data = df)
summary(two.way1)
# if density:fertilizer's Pr(>F)  < 0.05, then analyze with two.way1 instead with two.way
TukeyHSD(two.way)
#TukeyHSD(two.way1)
# Then, print with box-line picture
library(ggplot2)
ggplot(data=df, aes(x=fertilizer, y=yield)) +
    geom_boxplot() +
    facet_wrap(~density)
# aes(x= down, facet_wrap(~ up
# choose the hignest one


## 5. ggplot2 draw picture to show 单因素方差分析
df <- read.csv("example_data/anova_example.csv", header = T)
library(reshape2)
df1<-melt(df)
library(ggplot2)
ggplot(df1, aes(x=variable, y=value))+
    geom_boxplot(aes(fill=variable), show.legend = F)+
    theme_bw()+
    labs(x="aaa", y="")+
    annotate("text", x=1, y=630, label="a", size=5)+
    annotate("text", x=2, y=670, label="a", size=5)+
    annotate("text", x=3, y=710, label="a", size=5)+
    annotate("text", x=4, y=635, label="a", size=5)+
    scale_fill_manual(values = c("#7fd2ff", "#7fd2ff", "#ff9d1e", "#7fd2ff"))
