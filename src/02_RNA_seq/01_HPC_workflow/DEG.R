.libPaths("/public/home/ychang/anaconda3/lib/R/library")
library(DESeq2)
library(ggplot2)
library(magrittr)
library(dplyr)
library(openxlsx)
library(stringr)

##load data and contrast group

countdata <- read.csv("gene_count.csv", row.names="gene_id")
complist <- read.csv("contrast.csv", header=F,stringsAsFactors =F)
col <- read.csv("sample_info.csv", sep=",", row.names=1)
countdata <- countdata[,rownames(col)]
all(rownames(col) %in% colnames(countdata))
all(rownames(col) == colnames(countdata))
col$replicate <- as.factor(col$replicate)

##call overall DEG
DESeqDataSetFromMatrix(countData = countdata,colData = col, design =  ~replicate + condition) %>%
  DESeq(parallel=T) %>%
  vst() -> vst_all

#plot overall PCA
pca <- plotPCA(vst_all,intgroup=c("condition"),ntop=1000, returnData=T)
pdf("DEG_PCA.pdf", width =5,height =5,bg="transparent")
ggplot(pca,aes(PC1, PC2, color=condition)) +
  geom_point(size=3)+
  theme(axis.title=element_text(size=8,face = "bold"))+
  theme(axis.text=element_text(size=5,face = "bold",color = "black"))+
  theme(legend.text=element_text(size=4))+
  theme(legend.title =element_text(size=4))+
  xlab(paste0("PC1: ",round(100 * attr(pca, "percentVar"))[1],"% variance")) +
  ylab(paste0("PC2: ",round(100 * attr(pca, "percentVar"))[2],"% variance")) +
  geom_text(aes(label=name),hjust=0.5, vjust= 3, show.legend = F,cex=1.5,color="black")
dev.off()

#plot PCC heatmap
library(pheatmap)
pdf("DEG_PCC_heatmap.pdf", width =6,height =6,bg="transparent")
pheatmap(cor(as.matrix(assay(vst_all))),
         clustering_distance_rows=dist(cor(as.matrix(assay(vst_all)))),
         clustering_distance_cols=dist(cor(as.matrix(assay(vst_all)))),
         display_numbers=T,
         main="Pearson Correlation of DEGs between Samples",
         col=colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Oranges")) )(255))
dev.off()


#load anno list
anno <- read.xlsx("~/genome_anno/mh63/MH63RS3.encode.xlsx")

#call DEGs of desired groups separately

deg <- function(x){
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

  write.table(file=paste0("up_",x[1],"_vs_",x[2],".txt"),up,quote = F,sep="\t",row.names=F)
  write.table(file=paste0("down_",x[1],"_vs_",x[2],".txt"),down,quote = F,sep="\t",row.names=F)
  write.table(file=paste0("allsignif_",x[1],"_vs_",x[2],".txt"),all_deg,quote = F,sep="\t",row.names=F)
  write.table(file=paste0("allGenes_",x[1],"_vs_",x[2],".txt"),all,quote = F,sep="\t",row.names=F)
}
lapply(complist,deg)

save.image(file = "DEG.RData")