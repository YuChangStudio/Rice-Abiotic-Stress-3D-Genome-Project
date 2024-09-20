library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggnewscale)
library(openxlsx)
library(stringr)
library(ral)
library(pheatmap)

####===================================######
####Load overall configs
####===================================######

#color
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector
col.DS <- ralcolors["RAL8002"]%>%as.vector
col.RE <- ralcolors["RAL1014"]%>%as.vector

#ggplot2 theme
theme(panel.background = element_rect(fill = "transparent",colour='black',size = 1),
      panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
      panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
      plot.background = element_blank(),
      axis.ticks.length = unit(0.4,"lines"),
      axis.ticks = element_line(color='black',size =0.45 ),
      axis.line = element_line(colour = "black",size = 0),
      axis.title.x=element_text(colour='black', size=12,face = "bold"),
      axis.title.y=element_text(colour='black', size=12,face = "bold",angle = 90),
      axis.text.x=element_text(colour='black',size=10,face = "bold"),
      axis.text.y = element_text(color = "black",size = 10,face = "bold"),
      strip.text.x = element_text(colour='black',angle = 0,size = 10,face = "bold"),
      strip.text.y = element_text(colour='black',angle = -90,size = 10,face = "bold"),
      #legend.position=c(0.056,0.88),
      legend.title = element_text(colour='black',size = 8,face = "bold"),
      legend.background = element_rect(fill = "white",colour = "black"),
      legend.key = element_blank())%>%theme_set()

#Functional gene annotations derived from funRiceGenes (March, 2021).
anno <- read.xlsx("~/genome_anno/MH63RS3/MH63RS3_encode_Sep2022.xlsx")

####===================================######
####Step1: Data filtring and DEG calling by DEseq2
####===================================######

#Load data
countdata <- read.csv("../01_data_table_and_sample_list/gene_count.csv", row.names="gene_id")
complist <- read.csv("../01_data_table_and_sample_list/contrast.csv", header=F,stringsAsFactors =F)
col <- read.csv("../01_data_table_and_sample_list/sample_info.csv", sep=",", row.names=1)
countdata <- countdata[,rownames(col)]
all(rownames(col) %in% colnames(countdata))
all(rownames(col) == colnames(countdata))
col$replicate <- as.factor(col$replicate)

#Exclude genes with replicate-wise average TPM < 1 in all three conditions.
tpm <- read.table("../01_data_table_and_sample_list/gene.tpm",sep = "\t",header = T)%>%
  filter((MH_D_1+MH_D_2)/2 >=1 | (MH_N_1+MH_N_2)/2 >=1 | (MH_RE_1+MH_RE_2)/2 >=1)%>%arrange(gene_id)

countdata <- countdata[tpm$gene_id,]

#Calculate VST for QC plot
DESeqDataSetFromMatrix(countData = countdata,colData = col, design =  ~replicate + condition) %>%
  DESeq(parallel=T) %>%
  vst() -> vst_all

#Call DEGs of desired groups separately
call.deg <- function(x){
  col_sub <- col[col$condition==x[1]|col$condition==x[2],,drop=F]
  DESeqDataSetFromMatrix(countData = countdata[,rownames(col_sub)],
                         colData = col_sub,
                         design =  ~replicate + condition) %>%
    DESeq(parallel=T) %>%
    results(contrast = c("condition",x[1],x[2])) %>%
    as.data.frame() %>%
    merge(.,anno,by.x="row.names",by.y="MH63RS3_ID",all.x=T) %>%
    assign(paste0("res_",x[1],"_vs_",x[2]),.,envir=.GlobalEnv)-> all

  all_deg <- subset(all,padj<0.05)
  up <- subset(all,padj<0.05 & (log2FoldChange > 1))
  down <- subset(all,padj<0.05 & (log2FoldChange < -1))

  write.table(file=paste0("../02_DEG/01_DEG_table/up_",x[1],"_vs_",x[2],".txt"),up,quote = F,sep="\t",row.names=F)
  write.table(file=paste0("../02_DEG/01_DEG_table/down_",x[1],"_vs_",x[2],".txt"),down,quote = F,sep="\t",row.names=F)
  write.table(file=paste0("../02_DEG/01_DEG_table/allGenes_",x[1],"_vs_",x[2],".txt"),all,quote = F,sep="\t",row.names=F)

  write.xlsx(file = paste0("../02_DEG/02_DEG_spreadsheet/up_",x[1],"_vs_",x[2],".xlsx"),up)
  write.xlsx(file = paste0("../02_DEG/02_DEG_spreadsheet/down_",x[1],"_vs_",x[2],".xlsx"),down)
  write.xlsx(file = paste0("../02_DEG/02_DEG_spreadsheet/allGenes_",x[1],"_vs_",x[2],".xlsx"),all)

}
lapply(complist,call.deg)


####===================================######
####Step2: General QC
####===================================######
#plot VST-based Pearson correlation matrix
pdf("../02_DEG/03_DEG_plot/01_QC/VST_PCC.heatmap.pdf",height = 5,width = 5,bg=bg)
pheatmap(cor(as.matrix(assay(vst_all))),
         clustering_distance_rows=dist(cor(as.matrix(assay(vst_all)))),
         clustering_distance_cols=dist(cor(as.matrix(assay(vst_all)))),
         display_numbers=T,
         main="Pearson correlation in terms of VST",
         col=colorRampPalette( RColorBrewer::brewer.pal(9, "Oranges"))(255))
dev.off()

#plot VST-based PCA

pca <- plotPCA(vst_all,intgroup=c("condition"),ntop=1000, returnData=T)

pdf("../02_DEG/03_DEG_plot/01_QC/VST.pca.pdf",height = 2.5,width = 2.5,bg=bg)
ggplot(pca,aes(x=PC1,y=PC2))+
  geom_point(shape=21,color="black",aes(fill=group),size=4)+
  scale_fill_manual(values = c("MH_D"=col.DS,"MH_N"=col.NC,"MH_RE"=col.RE))+
  geom_text(aes(label=name))+
  xlab(paste0("Principal Component #1 [",round(100 * attr(pca, "percentVar"))[1],"%]")) +
  ylab(paste0("Principal Component #2 [",round(100 * attr(pca, "percentVar"))[2],"%]")) +
  theme(legend.position = "none")
dev.off()

####===================================######
####Step2: View DEGs by VocanoPlot
####===================================######

#plot DEG volcano plot.
#DS=drought stress; RE=recovery; NC=normal condition
deg.quant <- res_MH_D_vs_MH_N%>%select(log2FoldChange,padj)%>%mutate(cat="DS_vs_NC")
deg.quant <- res_MH_RE_vs_MH_D%>%select(log2FoldChange,padj)%>%mutate(cat="RE_vs_DS")%>%rbind(deg.quant,.)
deg.quant <- res_MH_RE_vs_MH_N%>%select(log2FoldChange,padj)%>%mutate(cat="RE_vs_NC")%>%rbind(deg.quant,.)

col.up <- ralcolors["RAL3014"]%>%as.vector
col.down <- ralcolors["RAL5012"]%>%as.vector
#col.ns <- ralcolors["RAL6021"]%>%as.vector

pdf(file = "../02_DEG/03_DEG_plot/MH_DEG.volcano.pdf",height = 2,width =6 ,bg = bg)
ggplot() +
  geom_hex(data = deg.quant%>%filter(abs(log2FoldChange) < 1 | padj >= 0.05),
           aes(x=log2FoldChange,y= -log10(padj)),bins=100)+
  scale_fill_gradient(low = "grey75" ,high = "black")+
  geom_point(data = deg.quant%>%filter(log2FoldChange >= 1 , padj < 0.05),
             aes(x=log2FoldChange,y= -log10(padj)),
             fill=alpha(col.up ,0.4),color=bg,shape=21,size=1,show.legend=c(fill=T,color=F,shape=F),inherit.aes = F) +
  geom_point(data = deg.quant%>%filter(log2FoldChange <= -1 , padj < 0.05),
             aes(x=log2FoldChange,y= -log10(padj)),
             fill=alpha(col.down,0.4),color=bg,shape=21,size=1,show.legend=c(fill=T,color=F,shape=F),inherit.aes = F)+
  facet_wrap(cat~.,nrow = 1,ncol = 3,scales = "free_x")
dev.off()

#plot data distribution in terms of log2FC
pdf(file = "../02_DEG/03_DEG_plot/MH_DEG.volcano.DataDistribute.pdf",height = 1.5,width =5 ,bg = bg)
ggplot() +
  geom_density(data = deg.quant%>%filter(log2FoldChange >= 1 , padj < 0.05),aes(x=log2FoldChange,y=after_stat(count)),fill=alpha(col.up ,0.5))  +
  geom_density(data = deg.quant%>%filter(log2FoldChange <= -1 , padj < 0.05),aes(x=log2FoldChange,y=after_stat(count)),fill=alpha(col.down ,0.5)) +
  facet_wrap(cat~.,nrow = 1,ncol = 3,scales = "free_x")
dev.off()

####===================================######
####Step3: Clustering genes
####===================================######
#Before plotting heatmap we perform hierarchy clustring seprately.
#Values are normalized by Z-transform
#Note the gene order in the hclust should be in consistent with that in the matrix for heatmap plotting

tpm.Z <- tpm %>% select(-gene_id)%>% +1%>%log2()%>%t()%>%scale()%>%t()%>%as.data.frame()

rownames(tpm.Z) <- tpm$gene_id

tpm.hclust <- dist(tpm.Z )%>%hclust()

#Cut the hclust result into 4 groups, annotate and export them.
tpm.hclust.cutree<- tpm.hclust %>% cutree(k = 4) %>% as.data.frame() %>% rename(h.cluster=".")

tpm.hclust.cutree <- tpm.hclust.cutree %>%merge(.,anno,by.x="row.names",by.y="MH63RS3_ID",all.x=T)

write.xlsx(file = "../02_DEG/04_hclust/hclust.geneInfo.xlsx",tpm.hclust.cutree,overwrite = T)


#Annotate the original TPM data table with DEG info

degUp.id.DS_vs_NC <- res_MH_D_vs_MH_N%>%filter(log2FoldChange >= 1 & padj < 0.05 )%>%select(Row.names)
degDown.id.DS_vs_NC <- res_MH_D_vs_MH_N%>%filter(log2FoldChange <= -1 & padj < 0.05 )%>%select(Row.names)

degUp.id.RE_vs_DS <- res_MH_RE_vs_MH_D%>%filter(log2FoldChange >= 1 & padj < 0.05 )%>%select(Row.names)
degDown.id.RE_vs_DS <- res_MH_RE_vs_MH_D%>%filter(log2FoldChange <= -1 & padj < 0.05 )%>%select(Row.names)

degUp.id.RE_vs_NC <- res_MH_RE_vs_MH_N%>%filter(log2FoldChange >= 1 & padj < 0.05 )%>%select(Row.names)
degDown.id.RE_vs_NC <- res_MH_RE_vs_MH_N%>%filter(log2FoldChange <= -1 & padj < 0.05 )%>%select(Row.names)

tpm <- tpm%>%mutate(DEG_DS_vs_NC=ifelse(gene_id %in% degUp.id.DS_vs_NC$Row.names,"up",
                                        ifelse(gene_id %in% degDown.id.DS_vs_NC$Row.names,"down","ns")))%>%
              mutate(DEG_RE_vs_DS=ifelse(gene_id %in% degUp.id.RE_vs_DS$Row.names,"up",
                                        ifelse(gene_id %in% degDown.id.RE_vs_DS$Row.names,"down","ns")))%>%
              mutate(DEG_RE_vs_NC=ifelse(gene_id %in% degUp.id.RE_vs_NC$Row.names,"up",
                                        ifelse(gene_id %in% degDown.id.RE_vs_NC$Row.names,"down","ns")))%>%
              arrange(gene_id)



####===================================######
####Step4: Use heatmap to illustrate gene expression profiles
####===================================######
library(viridis)

#row.names(tpm) <- tpm$gene_id
#tpm <- tpm%>%select(-gene_id) #Note that pheatmap requires the rownames to be set manually.

pdf(file = "../02_DEG/03_DEG_plot/MH_DEG.heatmap2.pdf",height = 7,width =5 ,bg = bg)
pheatmap(#tpm%>%select(-DEG_DS_vs_NC,-DEG_RE_vs_DS,-CEM)%>%+1%>%log2(),
         tpm.Z,
         col=inferno(n = 255,alpha = 1,direction = -1),
         cluster_rows = tpm.hclust,
         cutree_rows =4,show_rownames=F, #The heatmap is hierarchically clustered and then cut to 6 groups.
         cellwidth = 12,cellheight = 0.01,treeheight_row = 20,treeheight_col = 10,
         annotation_row=tpm%>%select(DEG_DS_vs_NC,DEG_RE_vs_DS),
         annotation_colors  = list(
                               #DEG_RE_vs_NC=c(ns=bg,up=col.up,down=col.down),
                               DEG_RE_vs_DS=c(ns=bg,up=col.up,down=col.down),
                               DEG_DS_vs_NC=c(ns=bg,up=col.up,down=col.down)
                               #CEM=c(black="black",blue="blue",brown="brown",green="green",
                                     #red="red",turquoise="turquoise",yellow="yellow",grey=bg)
                               ))
dev.off()

#####==========================================================================##########
##output DS dynamic genes (show opposite DE patterns under DS_vs_NC and RE_vs_DS) and all gene info.
#####==========================================================================##########

##Output the info for all genes
gene.cat <- read.table("../01_data_table_and_sample_list/gene_count.csv",sep = ",",header = T)%>%
  select(gene_id)%>%
  mutate(exped=ifelse(gene_id %in% tpm$gene_id,"Yes","No"))%>%
  mutate(DEG.DS_vs_NC=ifelse(gene_id %in% degUp.id.DS_vs_NC$Row.names,"Up",
                             ifelse(gene_id %in% degDown.id.DS_vs_NC$Row.names,"Down","NS")))%>%
  mutate(DEG.RE_vs_DS=ifelse(gene_id %in% degUp.id.RE_vs_DS$Row.names,"Up",
                             ifelse(gene_id %in% degDown.id.RE_vs_DS$Row.names,"Down","NS")))%>%
  mutate(Dynamic.to.DS=ifelse(DEG.DS_vs_NC=="Up" & DEG.RE_vs_DS=="Down","Up",
                              ifelse(DEG.DS_vs_NC=="Down" & DEG.RE_vs_DS=="Up","Down","None")))%>%
  rename(geneID=gene_id)%>%
  left_join(.,WGCNA.ME.info)

gene.cat$CEM[is.na(gene.cat$CEM)] <- "grey"

#gene.cat$CEM <- factor(gene.cat$CEM,levels = rev(c("black","green","turquoise","yellow","red","blue","brown","grey")),ordered = T)
gene.cat$CEM <- factor(gene.cat$CEM,levels = rev(c("turquoise","black","red","yellow","green","blue","brown","grey")),ordered = T)
gene.cat$Dynamic.to.DS <- factor(gene.cat$Dynamic.to.DS,levels = rev(c("Up","Down","None")),ordered = T)

write.xlsx(file = "../07_Dynamic_DEG_group/all_gene_categories.xlsx",gene.cat)


##Dynamic genes
DS.dyna.up.genes <- gene.cat%>%filter(Dynamic.to.DS=="Up")%>%left_join(.,anno%>%rename(geneID=MH63RS3_ID),by="geneID")%>%
  select(-exped, -DEG.DS_vs_NC, -DEG.RE_vs_DS, -Dynamic.to.DS)

DS.dyna.down.genes <- gene.cat%>%filter(Dynamic.to.DS=="Down")%>%left_join(.,anno%>%rename(geneID=MH63RS3_ID),by="geneID")%>%
  select(-exped, -DEG.DS_vs_NC, -DEG.RE_vs_DS, -Dynamic.to.DS)

write.xlsx(file = "../07_Dynamic_DEG_group/DS_dyna_up_genes.xlsx",DS.dyna.up.genes)
write.xlsx(file = "../07_Dynamic_DEG_group/DS_dyna_down_genes.xlsx",DS.dyna.down.genes)

write.xlsx(file = "../07_Dynamic_DEG_group/Expressed_gene_categories.xlsx",tpm)

