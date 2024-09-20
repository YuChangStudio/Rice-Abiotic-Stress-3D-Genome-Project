###General setting
##set working directory
setwd("/public/home/ychang/3D_project/01_ChIPseq/08_DBA/02_by_narrowPeak/V2/R_wd")

##load require package
#The script reuqires DiffBind version > 3.0
#the "ChIPseeker" pkg must be the github developer version, or "subset" will not work on the annotated peak object.
library(DiffBind)
library(magrittr)
library(GenomicFeatures)
library(ChIPseeker)
library(dplyr)
library(openxlsx)
library(stringr)
library(ggrepel)
library(ggplot2)
library(S4Vectors)

### Overall configs
## Set plot theme
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
        strip.text.x = element_text(angle = 0,size = 6,face = "bold"),
        #legend.position=c(0.056,0.88),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white",colour = "black"),
        legend.key = element_blank())%>%theme_set()

bg="transparent"

##Peak annotation configs
#options(ChIPseeker.ignore_1st_exon = TRUE)
#options(ChIPseeker.ignore_1st_intron = TRUE)
#options(ChIPseeker.ignore_downstream = TRUE)
#options(ChIPseeker.ignore_promoter_subcategory = TRUE)

## Load gene annotation
ricefungene <- read.xlsx("~/genome_anno/mh63/MH63RS3.encode.xlsx")%>%
                rename(geneId=MH63RS3_ID)%>%
                select(geneId,MSU7_ID,RAPdb_ID,Uniprot_ID,Gene_Symbol,Keyword,Title,Description)

## Load MH63 genomic feature
mh63rs3.txdb <- makeTxDbFromGFF(format = "gff3",file ="~/genome_anno/mh63/MH63RS3.gff3")



###DBA analysis

##Load samples
sample.list <- read.csv("../00_info/sample_info.csv")
samples.dba <- dba(sampleSheet = sample.list)

##Count reads within the peak regions (defined in sample_info.csv).
# Turn off summits option for H3K9ac occurs in a relative board region (summits=F).
# The input read counts are subtracted from the IP read counts for each site in each sample (bSubControl=T).
# If the input library size is larger than the IP, it will be sub-sampled to a matching size (bScaleControl=T).
# Filter regions with reads FPKM <1 (filter=1) to skip the quantification on large regions with few reads (likely false-positive regions).
# All raw read counts are maintained for use by dba.analyze, regardless of how the 'score' option is set.
mtx.count <- dba.count(samples.dba,bParallel=T,score=DBA_SCORE_RPKM_FOLD,summits=F,
                  bScaleControl=T,bSubControl=T,filter=1,minCount=0)

# record count info
count.info <- dba.show(mtx.count)
libsizes <- cbind(LibReads=count.info$Reads, FRiP=count.info$FRiP,PeakReads=round(count.info$Reads * count.info$FRiP))
rownames(libsizes) <- count.info$ID

write.xlsx(file= "../03_processQC/01_librarySize.xlsx",libsizes,overwrite = T)


## Normalization
# Use libsize normalization (the default method).
# Note, although RLE is native to DESeq2 analysis, any of these normalization methods can be used with either the DESeq2 or edgeR analysis,
# as DiffBind converts normalization factors to work appropriately in either DESeq2 or edgeR .
# Set background=T if RLE or RiP methods are to be used to avoid making assumption about the changes in binding patterns.
mtx.norm <- dba.normalize(mtx.count,normalize=DBA_NORM_LIB,method=DBA_DESEQ2,library=DBA_LIBSIZE_FULL,background=F)

# record normalization info
libsizes.norm <- cbind(FullLibSize=mtx.norm$norm$DESeq2$lib.sizes,
                       NormFacs=mtx.norm$norm$DESeq2$norm.facs,
                       NormLibSize=round(mtx.norm$norm$DESeq2$lib.sizes/mtx.norm$norm$DESeq2$norm.facs))
rownames(libsizes.norm) <- count.info$ID

write.xlsx(file= "../03_processQC/02_librarySize.normalized.xlsx",libsizes.norm,overwrite = T)

## Establishing a model design and contrast
# Contrast factor is set as DBA_CONDITION (Condition column) in 01_info_sample_info.csv
# Auto contrasts based on the factor is generated and analyzed for overall profiling

# In DiffBind3, the 'minMembers' must be specified manually, or automatic contrasts mode will not take effect.
mtx.contrast <- dba.contrast(mtx.norm,categories=DBA_CONDITION,minMembers=2)

## Performing the differential analysis on all combinations of groups specified in DBA_CONDITION.
mtx.analysis <- dba.analyze(mtx.contrast,method=DBA_DESEQ2,bParallel=T)

## Retrieving the differentially bound/enriched sites
mtx.report <- dba.report(mtx.analysis,method=DBA_DESEQ2,bUsePval=F,fold=1,th=0.05,bNormalized=T,bNotDB=T,bDB=T,bAll=T)

## Profiling based on the normalized read intensity in the peak regions.
#plot PCA
pdf(file="../04_plot/01_DBA_PCA.pdf",bg=bg,width = 7,height = 7)
dba.plotPCA(mtx.analysis, DBA_CONDITION,label=DBA_ID,report=mtx.report)
dev.off()

#plot PCC
PCC <- dba.plotHeatmap(DBA = mtx.analysis,distMethod="pearson",maxSites=5000,report=mtx.report)
pdf(file="../04_plot/02_DBA_PCC.heatmap.pdf",bg=bg,width = 10,height = 10)
dba.plotHeatmap(DBA = mtx.analysis,
                report=mtx.report,
                distMethod="pearson",
                maxSites=5000,
                cellnote=PCC,
                notecol="black",
                colScheme=colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges") )(255))
dev.off()

#

### The following analysis aims at desired comparing groups
## Get DBAs between desired groups.
complist <- read.csv(file="../00_info/contrast.csv",header=F,stringsAsFactors =F)

getDBA <- function(x){
  #list_sub <- slist[slist$Condition==x[1]|slist$Condition==x[2],]
  #samplelist_sub <- dba(sampleSheet = list_sub)
  #countmatrix <- dba.count(samplelist_sub,minOverlap=1, bParallel=T)

        dba.contrast(mtx.norm,design=F,group1=dba.mask(mtx.norm,DBA_CONDITION,x[1]),
                      group2=dba.mask(mtx.norm,DBA_CONDITION,x[2]),name1=x[1],name2=x[2]) %>%
        dba.analyze(method=DBA_DESEQ2,bParallel=T) %>%
        dba.report(method=DBA_DESEQ2,fold = 0,th = 1,bUsePval = F) %>%
        assign(paste0("report_",x[1],"_vs_",x[2]),.,envir=.GlobalEnv)%>%
        annotatePeak(TxDb = mh63rs3.txdb,tssRegion = c(-1000,100),flankDistance=1000,addFlankGeneInfo=T) %>%
        assign(paste0("anno_",x[1],"_vs_",x[2]),.,envir=.GlobalEnv)
}
lapply(complist,getDBA)

#output data and get DBA object
outAnnoRes <- function(x){
  res <- get(paste0(x))

  all.name <- str_replace(paste0(x),pattern = "anno_","dba.all.")
  up.name <- str_replace(paste0(x),pattern = "anno_","dba.up.")
  down.name <- str_replace(paste0(x),pattern = "anno_","dba.down.")

up<-subset(res, Fold > 1 & FDR < 0.05)%>%
    assign(paste0(up.name),.,envir = .GlobalEnv)%>%
    as.data.frame()%>%
    merge(.,ricefungene,by="geneId",all.x=T)

    write.table(file = paste0("../05_DBA_txt/",up.name,".txt"),up,sep = "\t",quote = F)
    write.xlsx(file = paste0("../06_DBA_xlsx/",up.name,".xlsx"),up,overwrite = T)

down<-subset(res, Fold < -1 & FDR < 0.05)%>%
    assign(paste0(down.name),.,envir = .GlobalEnv)%>%
    as.data.frame()%>%
    merge(.,ricefungene,by="geneId",all.x=T)

    write.table(file = paste0("../05_DBA_txt/",down.name,".txt"),down,sep = "\t",quote = F)
    write.xlsx(file = paste0("../06_DBA_xlsx/",down.name,".xlsx"),down,overwrite = T)

write.table(file = paste0("../05_DBA_txt/",all.name,".txt"),res,sep = "\t",quote = F)
write.xlsx(file = paste0("../06_DBA_xlsx/",all.name,".xlsx"),res,overwrite = T)
}

ls(pattern = "anno_")%>%lapply(., outAnnoRes)

#plot tar matrix
pro<- getPromoters(TxDb = mh63rs3.txdb,upstream = 3000,downstream = 3000,by = "gene")

tag.mtx.list <- list()

gettagmtxlist <- function(x){
  df.name <- paste0(x)
  df <- get(df.name )
  tag.mtx.list[[df.name]] <<- getTagMatrix(df@anno,windows = pro)
}

ls(pattern = "^dba.")%>%lapply(., gettagmtxlist)

size=ls(pattern = "^dba.")%>%length()

pdf(file = "../04_plot/03_DBA.profileplot.pdf",height = 1.5*size ,width = 4,bg=bg)
plotAvgProf(tag.mtx.list,
            xlim=c(-3000,3000),
            xlab = "Distance to TSS",
            ylab = "DBA Distribution Frequency",
            facet = "row",
            free_y=F) +
  theme(strip.text.y = element_text(size=4)) +
  scale_x_continuous(expand = c(0,0),position = "bottom",
                     breaks = seq(-3000,3000,by=1000))
  #coord_cartesian(ylim = c(0.0001,0.0007))
dev.off()

#plot feature propotions
plot.annobar <-function(x){
  df.name <- paste0(x)
  df <- get(df.name )
  plotAnnoBar(df)%>%
  ggsave(.,filename = paste0("../04_plot/03_Feature_proportion.",df.name,".annobar.pdf"),height = 4,width = 4,bg=bg,device = "pdf")
}

ls(pattern = "^dba.")%>%lapply(., plot.annobar)

save.image(file = "chip_dba.RData") #script test checkpoint
quit("no")

#plot DBA vocano
voc.in <- data.frame()
getVocIn <- function(x){
  df.name <- paste0(x)%>%str_replace(.,pattern = "anno_",replacement = "")
  df <- get(paste0(x))

voc.in <<-df@anno%>%as.data.frame()%>%
          #merge(.,ricefungene,by="geneId",all.x=T)%>%
          select(Fold,FDR)%>%
          mutate(cat=df.name,color=as.factor( ifelse( FDR < 0.05 & Fold > 1 ,"up",
                                                      ifelse(FDR < 0.05 & Fold < -1 ,"down","nc"))))%>%
          rbind(voc.in,.)
}
ls(pattern = "^anno")%>%lapply(.,getVocIn)

#voc.in$Gene_Symbol[is.na(voc.in$Gene_Symbol)]<-""

size=ls(pattern = "^anno")%>%length()
size=3*size

pdf(file = "02_plot/04_DBA.volcano.pdf", width = size,height =5,bg="transparent")
ggplot(voc.in,aes(Fold, -log10(FDR))) +
  geom_point(aes(fill=color),size=1.5,shape=21,color="transparent") +
  scale_fill_manual(values = alpha(c("down" = "#00bfff",
                                     "nc" = "grey75",
                                     "up"= "tomato"),0.4))+
  #scale_color_manual(values = c("n"="transparent","ed"="black"))+
  scale_y_continuous(expand = c(0,0.2))+
  #scale_x_continuous(expand = c(0,0))+
  #coord_cartesian(ylim = c(0,15))+
  facet_wrap(cat~.)+
  #geom_text_repel(data=voc.in[!voc.in$color =="nc",],
  #                aes(label=voc.in[!voc.in$color =="nc",]$Gene_Symbo,size=0.8),
  #                segment.size=0.2,show.legend = F,force=10)+
  theme(
        strip.text.x = element_text(angle = 0,size = 6,face = "bold")
        )
dev.off()

save.image(file = "chip_dba2.RData")

#End of the process
