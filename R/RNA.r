###
#* @Date: 2021-08-09 16:58:31
#* @LastEditors: Hwrn
#* @LastEditTime: 2021-08-09 17:25:24
#* @FilePath: /metaSC/R/RNA.r
#* @Description:
###

## First, for transcriptome, do mapping
# for sample in SRR11412240 SRR11412241 SRR11412242 SRR11412251 SRR11412252 SRR11412253
# do
#     in1=00_data/$sample/$sample.fastq
#     clean1=01_qc/${sample}.fq
#     map=02_map/${sample}

#     fastp \
#         -i ${in1} -o ${clean1} \
#         -h ${clean1}.fastp.html -j ${clean1}.fastp.json \
#         -q 20 -u 30 \
#         -w 20 \
#     >> Oerr/01_${sample}.log 2>&1

#     in1=00_data/$sample/$sample.fastq
#     clean1=01_qc/${sample}.fq
#     map=02_map/${sample}

#     hisat2 -x 00_data/grch37/genome \
#         -p 20 \
#         -U ${clean1} \
#         -S ${map}_raw.sam \
#     >> Oerr/02_${sample}.log 2>&1

#     samtools view -bS ${map}_raw.sam > ${map}_raw.bam
#     samtools sort ${map}_raw.bam -o ${map}_sort.bam \
#     >> Oerr/01_02_${sample}.log 2>&1

#     samtools index ${map}_sort.bam ${map}_sort.bam.bai \
#         -@ 20 \
#     >> Oerr/02_${sample}.log 2>&1
# done

# featureCounts -p \
#     -T 20 \
#     -t exon \
#     -g gene_id \
#     -a 00_data/gencode.v19.annotation.gff3.gz \
#     -o Analyze/counts.txt \
#     02_map/*_sort.bam \
#     >> Oerr/03.log 2>&1
### generate Analyze/counts.txt for further analysis

RawReadCounts <- read.table("Analyze/counts.txt",header = TRUE,row.names =1,quote="")
RawReadCounts <- as.matrix(RawReadCounts[,6:11])

library("DESeq2")
cl=as.factor(rep(c(0,1),c(3,3)))
coldata <- data.frame(row.names=colnames(RawReadCounts),cl)

dds <- DESeqDataSetFromMatrix(countData=RawReadCounts,colData=coldata,design=~cl)
(dds[ rowSums(counts(dds)) > 0, ])
dds1 <- dds[ rowSums(counts(dds)) > 1, ]
dds2 <- DESeq(dds1)
res_0.05 <- results(dds2,alpha=0.05)
res_0.05 <- res_0.05[order(res_0.05$padj),]
resdata_0.05 <- merge(as.data.frame(res_0.05),as.data.frame(counts(dds2,normalized=TRUE)),by="row.names",sort=F)

valueFind <- function(x, type = "gene_id"){
    for(i in 1:length(x)){
        if (startsWith(x[[i]], type)){
        return(strsplit(x[[i]], '=')[[1]][2])
        }
    }
    return(NA)
}

gff9=read.table('gencode.v19.annotation.gff3.gz', header = F, sep = '\t', stringsAsFactors = F)$V9
gff9 <- data.frame(do.call(rbind, strsplit(gff9, split = ";")), stringsAsFactors=F)
gff_name <- data.frame(gene_id = as.vector(apply(gff9,1,valueFind)),
                       gene_name = as.vector(apply(gff9,1,function(x){valueFind(x,type = "gene_name")})),
                       stringsAsFactors=F)
gff_name=gff_name[!duplicated(gff_name),]

resdata_0.05=resdata_0.05[!duplicated(resdata_0.05$Row.names),]
resdata_0.05=merge(gff_name, resdata_0.05, by.x='gene_id', by.y='Row.names',all.y=T,sort=T)
names(resdata_0.05)

(res.top100=resdata_0.05[order(abs(resdata_0.05$log2FoldChange), decreasing = T),][1:100,c(2,4)])

library(pheatmap)
for (signcutoff in c(1,1.5,2)){
    diff_gene_deseq2 <- subset(resdata_0.05, padj < 0.01 & abs(log2FoldChange) > signcutoff)
    png(paste(signcutoff, "_rnaseq.png",sep = ''), width=600, height=600)
    pheatmap(diff_gene_deseq2[9:14],main = paste('FoldChange cutoff: ', signcutoff),labels_row = diff_gene_deseq2$gene_name)
    dev.off()
}

## Micor-array
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)

library(affy)
affyData <- ReadAffy(celfile.path = "GSE17400/")
gn=geneNames(affyData)

deg <- AffyRNAdeg(affyData)
summaryAffyRNAdeg(deg)
plotAffyRNAdeg(deg)

library(simpleaffy)
Data.qc <- qc(affyData)
plot(Data.qc)

exprs(affyData)
CLLmas5 = mas5(affyData)
sampleNames(CLLmas5) <- gsub(".CEL.gz$", "", sampleNames(CLLmas5))

emat.mas5=as.data.frame(exprs(CLLmas5))
emat.mas5.log2=log2(emat.mas5)

results.mas5=data.frame((emat.mas5.log2[,c(1,4)]+emat.mas5.log2[,c(2,5)]+emat.mas5.log2[,c(3,6)])/3)
colnames(results.mas5)=c('mock', 'SARS')
results.mas5$log2FoldChange=results.mas5$SARS-results.mas5$mock

mas5.top100=results.mas5[order(abs(results.mas5$log2FoldChange), decreasing = T),][1:100,]
(mas5.top100=merge(ids, mas5.top100, by.x='probe_id', by.y=0, all.y=T, sort=T)[,c(2,5)])
overlap.top100=merge(res.top100, mas5.top100, by.x=1, by.y=1, all.x=T, all.y=T, sort=T)
write.csv(overlap.top100, file="Analyze/overlap.top100.csv", quote=F, row.names = F)

library(multtest)
library(siggenes)

cl=rep(c(0,1),c(3,3))
sam.out=sam(emat.mas5.log2,cl,gene.names = gn)

library(pheatmap)
for (signcutoff in c(1,1.5,2)){
    sum.sam.out=summary(sam.out,signcutoff)

    list.siggenes(sam.out,signcutoff)
    siggn=emat.mas5.log2[names(sum.sam.out@row.sig.genes),]

    siggn=merge(ids, siggn, by.x='probe_id', by.y=0, all.y=T, sort=T)
    png(paste(signcutoff, "_microarr.png",sep = ''), width=600, height=600)
    pheatmap(siggn[3:8],main = paste('signcutoff: ', signcutoff),labels_row = siggn$symbol)
    dev.off()
}
