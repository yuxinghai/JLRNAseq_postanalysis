setwd("/home/yuxh/hg19_jlRNAseq/")
directory <-"/home/yuxh/hg19_jlRNAseq/"
#distance in heatmap(all samples)
library(DESeq2)
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
#head(sampleFiles)
sampleCondition <-c("ipsf","ipsf","control","control","inono","inono","iNT","iNT","ipspc1"
                    ,"ipspc1","iV2","iV2")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition) 
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','ipsf',"inono","iNT"
                                                                          ,"ipspc1","iV2"))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-c("F1","F2","G1","G2","NO1","NO2","NT1","NT2","PS1","PS2","V2_1","V2_2")
res<-results(dds)
rld<- rlogTransformation(dds, blind=TRUE)
sampleDists <- dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:12, 1:12]
sampleDistMatrix <- as.matrix( sampleDists )
#colnames(sampleDistMatrix) <- NULL   
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
pdf("sample distance.pdf")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
Heatmap(sampleDistMatrix,cluster_rows =F,cluster_columns=F,col =colours ,name = "value",heatmap_legend_param = list(color_bar = "continuous"))
dev.off()
