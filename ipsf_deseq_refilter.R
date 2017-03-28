setwd("/home/yuxh/expres_fc/ipsf_control")

#F1,F2 on behalf of ipsf1 and ipsf2, G1 and G2 is iGFP on behalf of control
# gene different expression analysis
library(DESeq2)
#directory <- "/data3/zhoulab/yuxinghai/jlRNAseq/02_mapping/expres_data"
directory <-"/home/yuxh/expres_fc/ipsf_control"
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
#head(sampleFiles)
sampleCondition <-c("ipsf","ipsf","control","control")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','ipsf'))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-c("F1","F2","G1","G2")
# no filter the gene count number 
res<-results(dds)

#write all expression gene name to a file
all_exp <-res[res$baseMean !=0,]
write.table(row.names(all_exp),"ipsf_all_exp_name.tsv",sep="\n",quote = F,col.names = F,row.names = F)
res<-res[order(res$padj),]
mcols(res, use.names=TRUE)
head(res)

#add gene symbol and filter diff expressed gene
#up:log2FC>0.58, down:log2fc<-0.58
res <-res[complete.cases(res),] 
head(res)
res <-res[ res$padj <0.05,]
res$name <- rownames(res)
data_res <-as.data.frame(res)
dim(data_res)
genesymbol <- read.table("/home/yuxh/CSLab/mart_export.txt",header = T,sep = "\t")
#new is a dataframe with gene symbol
new <- merge(genesymbol,data_res,by.y="name",by.x="Ensembl.Gene.ID")
dim(new) #it's equal to dim(data_res)
down_gene <- new[new$log2FoldChange<(-0.58),]
#sort by -log2fc
down_gene <-down_gene[with(down_gene,order(-abs(log2FoldChange),padj)),]
up_gene <- new[new$log2FoldChange>log2(1.5),]
up_gene <-up_gene[with(up_gene,order(-abs(log2FoldChange),padj)),]
up_and_down<-new[(new$log2FoldChange<(-0.58)) | (new$log2FoldChange>log2(1.5)), ]
#genename for go term
write.table(up_and_down$Ensembl.Gene.ID,"ipsf_vs_control_sig_name.csv",sep = "\n",quote = F,col.names = F,row.names = F)
write.table(down_gene$Ensembl.Gene.ID,"ipsf_vs_control_down_name.csv",sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$Ensembl.Gene.ID,"ipsf_vs_control_up_name.csv",sep = "\n",quote = F,col.names = F,row.names = F)

write.csv(down_gene,"ipsf_vs_control_downgene.csv",quote = F,row.names = F)
write.csv(up_gene,"ipsf_vs_control_upgene.csv",quote = F,row.names = F)

#MAplot
plotMA(dds,ylim=c(-4,4),main='ipsf vs control')
dev.copy(png,'ipsf_vs_control_MAplot.png')
dev.off()

#heatmap
rld<- rlogTransformation(dds, blind=TRUE)

#association
pdf("correlation between F2 VS F1, G2 VS G1.pdf")
x1 <- assay(rld)[,1]
y1 <- assay(rld)[,2]
cor.test(x1,y1)

plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3,sub="R=0.99")
x2 <- assay(rld)[,3]
y2 <- assay(rld)[,4]
cor.test(x2,y2)
plot( assay(rld)[, 3:4], col="#00000020", pch=20, cex=0.3 ,sub="R=0.99")
#F1 and G2, F2 and G1
cor.test(x1,y2)
cor.test(x2,y1)
plot( assay(rld)[, c(1,4)], col="#00000020", pch=20, cex=0.3 ,sub="R=0.96")
plot( assay(rld)[, 2:3], col="#00000020", pch=20, cex=0.3 ,sub="R=0.96")
dev.off()


#all completed data
pdf('ipsf_vs_control_heatmap.pdf')
par(cex.lab=0.3)
res<-results(dds)
data <-res[ res$baseMean !=0,]
data <-data[complete.cases(data),] 
data<-data[order(-data$log2FoldChange),]
data$name <-row.names(data)
#change DESeqResults to dataframe
data <-as.data.frame(data)


raw_count <-as.data.frame(assay(dds))
raw_count <-cbind(name=row.names(assay(dds)),raw_count)
#gene_count is a merge data include rawcount and deseq2result
gene_count <-merge(data,raw_count,by="name")
gene_count_s <-gene_count[with(gene_count,order(-log2FoldChange)),]
# different expression data,get upgene and downgene from former text
Sig_change<-gene_count_s$name %in% up_and_down$Ensembl.Gene.ID



notSig_Change <-gene_count_s[!Sig_change,]
notSig_Change_Sort <-notSig_Change[with(notSig_Change,order(log2FoldChange)),][8:11]
dim(notSig_Change_Sort)
up_bool <-gene_count_s$name %in% up_gene$Ensembl.Gene.ID
up <-gene_count_s[up_bool,]
up <-up[with(up,order(-log2FoldChange)),][8:11]
dim(up)
down_bool<-gene_count_s$name %in% down_gene$Ensembl.Gene.ID
down <-gene_count_s[down_bool,]
down <-down[with(down,order(-log2FoldChange)),][8:11]
dim(down)
library("edgeR")
#draw heatmap order by -logfc
pre_draw <-rbind(up,notSig_Change_Sort,down)
pre_draw <-cpm(pre_draw)
pre_draw_scaled = t(apply(pre_draw, 1, scale))

library("gplots")
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
split = data.frame(c(rep("up",1366),rep("not_sig",16460),rep("down",1174)),
                   c(rep("1366",1366),rep("16460",16460),rep("1174",1174)))

ha = HeatmapAnnotation(df = data.frame(type=colnames(pre_draw)))
Heatmap(pre_draw_scaled,cluster_rows =F,show_row_names=F,split = split,
        col =colorRamp2(c(-0.8, 0, 0.8), c("green", "black", "red")),name = "value", 
        top_annotation = ha,heatmap_legend_param = list(color_bar = "continuous"),
        combined_name_fun = function(x) paste(x, collapse = "\n"))

dev.off()


#go term analysis

#library(biomaRt)
#genes = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl",version = 84)
# 得到go信息和gene
#gene2goInfo <- getBM(attributes=c('ensembl_gene_id','go_id','entrezgene','name_1006','go_linkage_type','namespace_1003'), mart = genes)
#write.table(gene2goInfo,"GENE_GO_ENSEMBLE84.tsv",quote = F,sep = "\t",row.names = F)
gene2goInfo <-read.table("GENE_GO_ENSEMBLE84.tsv",header = T,sep="\t")
# 过滤
# my own download from ensemble
gene2goInfo=gene2goInfo[gene2goInfo$go_id != "", ]
# 上述是一个gene对应一个go id，需要合并为一个gene对应多个go id，利用by函数（神器）
geneID2GO = by(gene2goInfo$go_id, gene2goInfo$ensembl_gene_id, function(x) as.character(x))
gene_go <-read.table("/home/yuxh/GO/ensemble84_gene_go.txt",header = T,sep = "\t")
go_2<-gene_go[gene_go$GO.Term.Accession !="",]
geneID2GO_2= by(gene_go$GO.Term.Accession, gene_go$Ensembl.Gene.ID, function(x) as.character(x))
# 差异表达基因
interesting_genes=factor(up_and_down$Ensembl.Gene.ID)
# 所有基因
all_genes <- data$name
# 构建基因列表
geneList <- factor(as.integer (all_genes %in% interesting_genes))
names(geneList)=all_genes
library("topGO")
library("org.Hs.eg.db")
library("GO.db")

MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
MF_resultFisher <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher")
MF_allRes <- GenTable(MF_GOdata,elimFisher = MF_resultFisher, orderBy = "elimFisher", topNodes = 10)

BP_GOdata<- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
BP_resultFisher <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher")
BP_allRes <- GenTable(BP_GOdata,elimFisher = BP_resultFisher, orderBy = "elimFisher", topNodes = 10)
CC_GOdata<- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
CC_resultFisher <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher")
CC_allRes <- GenTable(CC_GOdata,elimFisher = CC_resultFisher, orderBy = "elimFisher", topNodes = 10)
go_data <-rbind(MF_allRes,BP_allRes,CC_allRes)
go_data$go_type<-c(rep("MF",10),rep("BP",10),rep("CC",10))
library(ggplot2)
library(stringr)
go <-ggplot(go_data)
type=factor(go_data$go_type)
Pvalue=as.numeric(go_data$elimFisher)
go_term <-c(letters,"za","zb","zc","zd")
P <-go+geom_bar(aes(x=go_term,y=-log10(Pvalue),fill=go_type),stat="identity")+
  scale_x_discrete(labels =function(x) str_wrap(go_data$Term, width = 80))+coord_flip()+theme(axis.text.y =element_text(size=10))
gene_go <-read.table("/home/yuxh/GO/ensemble84_gene_go.txt",header = T,sep = "\t")
go_2<-gene_go[gene_go$GO.Term.Accession !="",]
geneID2GO_2= by(go_2$GO.Term.Accession, go_2$Ensembl.Gene.ID, function(x) as.character(x))
# 差异表达基因
interesting_genes=factor(sig$name)
# 所有基因
all_genes <- data$name
# 构建基因列表
geneList <- factor(as.integer (all_genes %in% interesting_genes))
names(geneList)=all_genes
MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
MF_resultFisher <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher")
MF_allRes <- GenTable(MF_GOdata,elimFisher = MF_resultFisher, orderBy = "elimFisher", topNodes = 10)

BP_GOdata<- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
BP_resultFisher <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher")
BP_allRes <- GenTable(BP_GOdata,elimFisher = BP_resultFisher, orderBy = "elimFisher", topNodes = 10)
CC_GOdata<- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
CC_resultFisher <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher")
CC_allRes <- GenTable(CC_GOdata,elimFisher = CC_resultFisher, orderBy = "elimFisher", topNodes = 10)
go_data <-rbind(MF_allRes,BP_allRes,CC_allRes)
#original term with imcomplete name

go_data$go_type<-c(rep("MF",10),rep("BP",10),rep("CC",10))
library(ggplot2)
library(stringr)
go <-ggplot(go_data)
type=factor(go_data$go_type)
Pvalue=as.numeric(go_data$elimFisher)
go_term <-c(letters,"za","zb","zc","zd")
P <-go+geom_bar(aes(x=go_term,y=-log10(Pvalue),fill=go_type),stat="identity")+
  scale_x_discrete(labels =go_data$Term)+coord_flip()




