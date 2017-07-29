
#go term analysis
#library(biomaRt)
#genes = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl",version = 74)
# 得到go信息和gene
#gene2goInfo <- getBM(attributes=c('ensembl_gene_id','go_id','entrezgene','name_1006','go_linkage_type','namespace_1003'), mart = genes)
#write.table(gene2goInfo,"/home/yuxh/GO/GENE_GO_ENSEMBLE74.tsv",quote = F,sep = "\t",row.names = F)

top_go_plot <-function(i_gene,pdf_name,expresed_gene=expresed_gene) {
  library("topGO")
  library("org.Hs.eg.db")
  library("GO.db")
  library(stringr)
  gene_go <-read.table("/home/yuxh/GO/GENE_GO_ENSEMBLE74.tsv",header = T,sep = "\t")
  #head(gene_go)
  go_2<-gene_go[gene_go$ensembl_gene_id !="",]
  geneID2GO_2= by(go_2$go_id, go_2$ensembl_gene_id, function(x) as.character(x))
  # 差异表达基因
  gene_id<-as.vector(i_gene$Ensembl.Gene.ID)
  intere_gene<-str_extract(gene_id,"ENSG[0-9]+")
  
  interesting_genes=factor(intere_gene)
  # 所有基因
  all_genes <- factor(str_extract(expresed_gene$name,"ENSG[0-9]+"))
  # 构建基因列表
  geneList <- factor(as.integer (all_genes %in% interesting_genes))
  names(geneList)=all_genes
  MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
  MF_resultFisher <- runTest(MF_GOdata, algorithm = "elim", statistic = "fisher")
  MF_allRes <- GenTable(MF_GOdata,elimFisher = MF_resultFisher, orderBy = "elimFisher", topNodes = 10,numChar=100)
  
  BP_GOdata<- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
  BP_resultFisher <- runTest(BP_GOdata, algorithm = "elim", statistic = "fisher")
  BP_allRes <- GenTable(BP_GOdata,elimFisher = BP_resultFisher, orderBy = "elimFisher", topNodes = 10,numChar=100)
  CC_GOdata<- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO_2)
  CC_resultFisher <- runTest(CC_GOdata, algorithm = "elim", statistic = "fisher")
  CC_allRes <- GenTable(CC_GOdata,elimFisher = CC_resultFisher, orderBy = "elimFisher", topNodes = 10,numChar=100)
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
    scale_x_discrete(labels =go_data$Term)+theme(axis.text.x = element_text(size =10,angle = 45))+coord_flip()
  P
  ggsave(paste(pdf_name,"_go.pdf",collapse=""),width = 15,height = 8.5)
  
}


de_and_qc <- function(directory,treat) {
  #head(sampleFiles)
  sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
  
  if (treat > "G1") {
    
    sampleCondition <-c("control","control",treat,treat)
    lev=c('control',treat)
  } else{
    
    sampleCondition <-c(treat,treat,"control","control")
    lev=c(treat,'control')
  }
  library(DESeq2)
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design= ~ condition)
  colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=lev)
  
  dds<-DESeq(ddsHTSeq)
  if (treat > "G1") {
    colnames(dds) <-c("G1","G2",paste(treat,1,collapse = ""),paste(treat,2,collapse = ""))
    
  } else{
    colnames(dds) <-c(paste(treat,1,collapse = ""),paste(treat,2,collapse = ""),"G1","G2")
  }
  
  
  
  # no filter the gene count number 
  res<-results(dds)
  write.table(row.names(res),paste(treat,"_all_name.tsv",collapse = ""),sep="\n",quote = F,col.names = F,row.names = F)
  
  
  #write all expression gene name to a file
  all_exp <-res[res$baseMean !=0,]
  write.table(row.names(all_exp),paste(treat,"_expresed_name.tsv",collapse = ""),sep="\n",quote = F,col.names = F,row.names = F)
  res<-res[order(res$padj),]
  
  #add gene symbol and filter diff expressed gene
  #up:log2FC>0.58, down:log2fc<-0.58
  res <-res[complete.cases(res),] 
  head(res)
  res <-res[ res$padj <0.05,]
  res$name <- rownames(res)
  data_res <-as.data.frame(res)
  dim(data_res)#[1] 7493    7
  
  genesymbol <- read.table("/home/yuxh/hg19_jlRNAseq/F/gencode.v19.annotation.gene_name",header = F,sep = " ")
  #new is a dataframe with gene symbol
  colnames(genesymbol) <-c("Ensembl.Gene.ID","genename")
  new <- merge(genesymbol,data_res,by.y="name",by.x="Ensembl.Gene.ID")
  dim(new) #8231    8
  
  down_gene <- new[new$log2FoldChange<(-0.58),]
  #sort by -log2fc
  down_gene <-down_gene[with(down_gene,order(-abs(log2FoldChange),padj)),]
  up_gene <- new[new$log2FoldChange>log2(1.5),]
  up_gene <-up_gene[with(up_gene,order(-abs(log2FoldChange),padj)),]
  up_and_down<-new[(new$log2FoldChange<(-0.58)) | (new$log2FoldChange>log2(1.5)), ]
  #genename for go term
  up_title=paste(treat,"_vs_control_up_gene.csv",collapse = "")
  down_title=paste(treat,"_vs_control_down_gene.csv",collapse = "")
  write.csv(down_gene,down_title,quote = F,row.names = F)
  write.csv(up_gene,up_title,quote = F,row.names = F)
  
  
  
  #MAplot
  plotMA(dds,ylim=c(-4,4),main=paste(treat,"_vs_control",collapse = ""))
  dev.copy(png,paste(treat,"_vs_control_MAplot.png",collapse = ""))
  dev.off()
  
  #heatmap
  rld<- rlogTransformation(dds, blind=TRUE)
  
  ###################association####################
  pdf(paste(treat,"correlation between G2 VS G1.pdf",collapse = ""))
  x1 <- assay(dds)[,1]
  y1 <- assay(dds)[,2]
  cor.test(x1,y1)
  #ata:  x1 and y1
  #t = 22630, df = 57818, p-value < 2.2e-16
  #alternative hypothesis: true correlation is not equal to 0
  #95 percent confidence interval:
  # 0.9999426 0.9999445
  #sample estimates:
  #  cor 0.9999436 
  #
  
  
  plot( assay(dds)[, 1:2], col="#00000020", pch=20, cex=0.3)
  x2 <- assay(dds)[,3]
  y2 <- assay(dds)[,4]
  cor.test(x2,y2)
  plot( assay(dds)[, 3:4], col="#00000020", pch=20, cex=0.3)
  #F1 and G2, F2 and G1
  cor.test(x1,y2)
  cor.test(x2,y1)
  plot( assay(dds)[, c(1,4)], col="#00000020", pch=20, cex=0.3 ,sub="R=0.96")
  plot( assay(dds)[, 2:3], col="#00000020", pch=20, cex=0.3 ,sub="R=0.96")
  dev.off()
  
  
  #all completed data
  
  res<-results(dds)
  expresed_gene <-res[ res$baseMean !=0,]
  expresed_gene <-expresed_gene[complete.cases(expresed_gene),] 
  expresed_gene<-expresed_gene[order(-expresed_gene$log2FoldChange),]
  expresed_gene$name <-row.names(expresed_gene)
  #change DESeqResults to dataframe
  expresed_gene <-as.data.frame(expresed_gene)
  
  
  raw_count <-as.data.frame(assay(dds))
  raw_count <-cbind(name=row.names(assay(dds)),raw_count)
  #gene_count is a merge data include rawcount and deseq2result
  gene_count <-merge(expresed_gene,raw_count,by="name")
  gene_count_s <-gene_count[with(gene_count,order(-log2FoldChange)),]
  # different expression data,get upgene and downgene from former text
  Sig_change<-gene_count_s$name %in% up_and_down$Ensembl.Gene.ID
  
  
  
  notSig_Change <-gene_count_s[!Sig_change,]
  notSig_Change_Sort <-notSig_Change[with(notSig_Change,order(log2FoldChange)),][8:11]
  dim(notSig_Change_Sort)
  up_bool <-gene_count_s$name %in% up_gene$Ensembl.Gene.ID
  up <-gene_count_s[up_bool,]
  up <-up[with(up,order(-log2FoldChange)),][8:11]
  up_num<-dim(up)[1]
  
  down_bool<-gene_count_s$name %in% down_gene$Ensembl.Gene.ID
  down <-gene_count_s[down_bool,]
  down <-down[with(down,order(-log2FoldChange)),][8:11]
  down_num<-dim(down)[1]
  
  
  ######go term###############
  top_go_plot(up_and_down,"up_and_Down",expresed_gene)
  top_go_plot(down_gene,"down",expresed_gene)
  top_go_plot(up_gene,"up",expresed_gene)
  
 
  library("edgeR")
  #draw heatmap order by -logfc
  pre_draw <-rbind(up,notSig_Change_Sort,down)
  pre_draw <-cpm(pre_draw)
  pre_draw_scaled = t(apply(pre_draw, 1, scale))
  all_num<-dim(pre_draw_scaled)[1]
  rest_num<-all_num-up_num-down_num

  library("gplots")
  library(ComplexHeatmap)
  library(circlize)
  library("RColorBrewer")
  pdf(paste(treat,'_vs_control_heatmap.pdf',collapse = ""))
  par(cex.lab=0.3)
  split = data.frame(c(rep("up",up_num),rep("not_sig",rest_num),rep("down",down_num)),
                     c(rep(as.character(up_num),up_num),rep(as.character(rest_num),rest_num),rep(as.character(down_num),down_num)))
  
  ha = HeatmapAnnotation(df = data.frame(type=colnames(pre_draw)))
  Heatmap(pre_draw_scaled,cluster_rows =F,show_row_names=F,split = split,
          col =colorRamp2(c(-0.8, 0, 0.8), c("green", "black", "red")),name = "value", 
          top_annotation = ha,heatmap_legend_param = list(color_bar = "continuous"),
          combined_name_fun = function(x) paste(x, collapse = "\n"))
  #dev.off in function will can't make the plot
 
  
}


########################################################################


