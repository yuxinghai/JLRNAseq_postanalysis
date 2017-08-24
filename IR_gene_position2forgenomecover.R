#############  plot  #############
setwd("~/PSF-INTON-RETATION/IRfinder/gene_position/")
prepare_for_plot<-function(filename){
  name2<-read.table(paste("~/PSF-INTON-RETATION/IRfinder/gene_position/",
                          filename,"_pre_mygenome.txt",sep=""),
                        header=F,stringsAsFactors = F,comment.char="#")[,c(4,5)]
  colnames(name2)<-c("position","height")
  name2$ratio <-name2$height/sum(name2$height)
  name2$ratio
}
All<-prepare_for_plot("KD_SS_PSF-v-Control.tab")
L10K<-prepare_for_plot("psf_IR_L10k")
L20K<-prepare_for_plot("psf_IR_L20k")
L30K<-prepare_for_plot("psf_IR_L30k")
L40K<-prepare_for_plot("psf_IR_L40k")
L50K<-prepare_for_plot("psf_IR_L50k")
L60K<-prepare_for_plot("psf_IR_L60k")
position <-seq(1,100)
comb<-cbind(All,L10K,L20K,L30K,L40K,L50K,L60K)

library(ggplot2)
library(reshape2)
m_comb <-melt(comb)
head(m_comb)
colnames(m_comb)<-c("position","type","ratio")
m_comb$type
ggplot(m_comb)+geom_line(aes(x=position,y=ratio,group=type ,col=type))+
  labs(y="percentage",title="Distribution of reteined intron")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("red", "yellow", "grey","green","blue","black","purple","pink"))
ggsave("~/PSF-INTON-RETATION/IRfinder/gene_position/IR_VS_control.intron_distribution.pdf")
