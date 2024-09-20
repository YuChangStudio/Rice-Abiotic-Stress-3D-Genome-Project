library(dplyr)
library(ggplot2)
library(pheatmap)
library(ral)
library(stringr)
#load plot theme
theme(panel.background = element_rect(fill = "transparent",colour='black',size = 1),
      panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
      panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
      axis.ticks.length = unit(0.4,"lines"),
      axis.ticks = element_line(color='black',size =0.45 ),
      axis.line = element_line(colour = "black",size = 0),
      axis.title.x=element_text(colour='black', size=12,face = "bold"),
      axis.title.y=element_text(colour='black', size=12,face = "bold",angle = 90),
      axis.text.x=element_text(colour='black',size=8),
      axis.text.y = element_text(color = "black",size = 8),
      strip.text.x = element_text(angle = 0,size = 12,face = "bold"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank())%>%theme_set()

bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector
col.DS <- ralcolors["RAL8002"]%>%as.vector
col.RE <- ralcolors["RAL1014"]%>%as.vector

condition.col<- c("NC"=col.NC,"DS"=col.DS,"RE"=col.RE)

###===============================================================================###
###calculate the PCC between replicates based on scaled read counts in 10kb bins.
###===============================================================================###
read.count.raw <- read.table("../01_bam_count/10kbp_bin/ChIApet.bam.count",header = T)

scale.factor <- read.table("../01_bam_count/10kbp_bin/ChIApet.bam.scalingFactors",header = T,row.names = 1)%>%
                t()%>%as.data.frame()

read.count.scale <- read.count.raw %>%
                    mutate(MH_H3K9ac_D_rep1.sorted.bam=MH_H3K9ac_D_rep1.sorted.bam * scale.factor$MH_H3K9ac_D_rep1.sorted.bam,
                           MH_H3K9ac_D_rep2.sorted.bam=MH_H3K9ac_D_rep2.sorted.bam * scale.factor$MH_H3K9ac_D_rep2.sorted.bam,
                           MH_H3K9ac_N_rep1.sorted.bam=MH_H3K9ac_N_rep1.sorted.bam * scale.factor$MH_H3K9ac_N_rep1.sorted.bam,
                           MH_H3K9ac_N_rep2.sorted.bam=MH_H3K9ac_N_rep2.sorted.bam * scale.factor$MH_H3K9ac_N_rep2.sorted.bam,
                           MH_H3K9ac_RE_rep1.sorted.bam=MH_H3K9ac_RE_rep1.sorted.bam * scale.factor$MH_H3K9ac_RE_rep1.sorted.bam,
                           MH_H3K9ac_RE_rep2.sorted.bam=MH_H3K9ac_RE_rep2.sorted.bam * scale.factor$MH_H3K9ac_RE_rep2.sorted.bam)

##plot PCC heatmap.
pdf("../03_plot/01_MH_H3K9ac_ChIApet.normReads.10kb.PCC.heatmap.pdf",height = 8,width = 8,bg=bg)
cor(read.count.scale%>%select(-chr,-start,-end),method = "pearson")%>%pheatmap(display_numbers = T)
dev.off()

##plot replicate-wise PCC scatter plot.
pcc.scatter.in <- rbind(read.count.scale%>%select(MH_H3K9ac_D_rep1.sorted.bam,MH_H3K9ac_D_rep2.sorted.bam)%>%
                          mutate(condition="DS")%>%rename(rep1=MH_H3K9ac_D_rep1.sorted.bam,rep2=MH_H3K9ac_D_rep2.sorted.bam),
                        read.count.scale%>%select(MH_H3K9ac_N_rep1.sorted.bam,MH_H3K9ac_N_rep2.sorted.bam)%>%
                          mutate(condition="NR")%>%rename(rep1=MH_H3K9ac_N_rep1.sorted.bam,rep2=MH_H3K9ac_N_rep2.sorted.bam),
                        read.count.scale%>%select(MH_H3K9ac_RE_rep1.sorted.bam,MH_H3K9ac_RE_rep2.sorted.bam)%>%
                          mutate(condition="RE")%>%rename(rep1=MH_H3K9ac_RE_rep1.sorted.bam,rep2=MH_H3K9ac_RE_rep2.sorted.bam))


pdf("../03_plot/PCC.pdf",height = 9,width = 3,bg=bg)
ggplot(pcc.scatter.in,aes(x=log2(rep1),y=log2(rep2)))+
  geom_point(aes(color=condition),size=0.5)+
  scale_color_manual(values = alpha(conditon.col,0.4))+
  geom_smooth(method = "lm",color="black",se=F)+
  xlab("Normalized reads abundance in 1kb bins of Repliact 1 (log2)")+
  ylab("Normalized reads abundance in 1kb bins of Repliact 2 (log2)")+
  facet_wrap(facets = .~condition,ncol = 1,scales = "free")+
  theme(legend.position = "None")
dev.off()



chiapet.prin=princomp(read.count.scale%>%select(-chr,-start,-end),cor=T)
summary(chiapet.prin)
chiapet.comp1=chiapet.prin$loadings[,1]
chiapet.comp2=chiapet.prin$loadings[,2]

chiapet.prin.pca<-cbind(as.data.frame(chiapet.comp1),as.data.frame(chiapet.comp2))

chiapet.prin.pca$sample <- rownames(chiapet.prin.pca)

library(ggforce)
pdf(file = "../03_plot/03_ChIA-PET.PCA.pdf",height = 3,width = 6,bg = bg)
ggplot(chiapet.prin.pca,aes(x=chiapet.comp1,y=chiapet.comp2,col=sample))+
  geom_point()+
  labs(x="PC1 (90.04%)",y="PC2 (6.55%)")
  #geom_rug()+
  #theme_bw()+
  #geom_mark_ellipse(aes(label=sample))
dev.off()