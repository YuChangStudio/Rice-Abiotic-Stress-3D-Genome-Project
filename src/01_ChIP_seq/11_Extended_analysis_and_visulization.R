library(DiffBind)
library(GenomicFeatures)
library(ChIPseeker)
library(openxlsx)
library(stringr)
library(ggrepel)
library(ggplot2)
library(S4Vectors)
library(hexbin)
library(ggnewscale)
library(ral)
library(ggsignif)
library(dplyr)

####=========================================================================================================######
####Load overall configs and proto-data done on HPC
####=========================================================================================================######

#ggplot theme
theme(aspect.ratio=1/1,

      line=element_line(size = 0.2351852, lineend = "square",color = "black"),
      rect=element_rect(fill=NULL,colour = "black",size = 0.2351852), #0.2351852mm=0.5pt under MACOS.
      text = element_text(family = "ArialMT",colour = "black"),

      panel.background = element_rect(fill = NULL,colour = "black"),
      panel.border = element_blank(),

      #panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
      #panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
      panel.grid = element_blank(), #journals such as plant cell do not like this...
      #axis.line = element_line(colour = "black",size = 0.5),
      axis.line = element_blank(),
      axis.ticks.length = unit(0.25,"lines"),
      axis.title.x=element_text(size=8),
      axis.title.y=element_text(size=8,angle = 90),
      axis.text=element_text(size=6),
      strip.text= element_text(size = 8),
      legend.title = element_blank(),
      legend.key = element_blank())%>%theme_set()

#color
bg="transparent"

#show_ralcolors() #show ral color map in IDE

col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1014"]%>%as.vector #"#DFCEA1"

condition.col<- c("NC"=col.NC,"DS"=col.DS,"RE"=col.RE)

load("chip_dba.RData") #proto data made on HPC through the dba3.R script (DiffBind version: 3.4.11).

ricefungene <- read.xlsx("~/genome_anno/MH63RS3/MH63RS3_encode_Sep2022.xlsx")%>%rename(geneId = MH63RS3_ID)

mh63rs3.txdb <- makeTxDbFromGFF(format = "gff3",file ="~/genome_anno/MH63RS3/MH63RS3.gff3")


####=========================================================================================================######
####Biological replicates QC based on the Pearson correlations of 10k-bined reads
####=========================================================================================================######

#The bined reads were count by muiltbamsummary from the deeptools toolkit.

#Normalize the read counts based on the scaling factors

scaleFactor <- read.delim("../09_BIned_read_count/bam.10k.scalingFactors",header = T,sep = "\t",row.names = 1)%>%
                t()%>%as.data.frame()

readCount.10kb <- read.delim("../09_BIned_read_count/bam.10k.count",header = T,sep = "\t")%>%
                  mutate(MH_D_H3K9ac_1_IP.rmdup.bam=MH_D_H3K9ac_1_IP.rmdup.bam*scaleFactor$MH_D_H3K9ac_1_IP.rmdup.bam,
                         MH_D_H3K9ac_2_IP.rmdup.bam=MH_D_H3K9ac_2_IP.rmdup.bam*scaleFactor$MH_D_H3K9ac_2_IP.rmdup.bam,

                         MH_N_H3K9ac_1_IP.rmdup.bam=MH_N_H3K9ac_1_IP.rmdup.bam*scaleFactor$MH_N_H3K9ac_1_IP.rmdup.bam,
                         MH_N_H3K9ac_2_IP.rmdup.bam=MH_N_H3K9ac_2_IP.rmdup.bam*scaleFactor$MH_N_H3K9ac_2_IP.rmdup.bam,

                         MH_RE_H3K9ac_1_IP.rmdup.bam=MH_RE_H3K9ac_1_IP.rmdup.bam*scaleFactor$MH_RE_H3K9ac_1_IP.rmdup.bam,
                         MH_RE_H3K9ac_2_IP.rmdup.bam=MH_RE_H3K9ac_2_IP.rmdup.bam*scaleFactor$MH_RE_H3K9ac_2_IP.rmdup.bam,

                         MH_input.rmdup.bam=MH_input.rmdup.bam*scaleFactor$MH_input.rmdup.bam)%>%
                  arrange(chr,start)

cor(readCount.10kb%>%select(-chr,-start,-end))

#The general PPCs between replicates are > 0.99.

####=========================================================================================================######
####Basic analysis based on the IDR peaks under each treatment.
####=========================================================================================================######

#profile peak width under DS (Drought Stress), NC (Normal Condition), and RE (REcovery).
#Use the non-redundent IDR peaks as input.

IDRpeaks.H3K9ac.NC <- read.table("../01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_N_H3K9ac.IDR.uniq.narrowPeak",sep = "\t")%>%
  select(V1,V2,V3)%>%mutate(peak.id=str_c("NC_peak_",row.names(.)),width=V3-V2,treatment="NC")%>%unique()%>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,seqnames.field = "V1",start.field = "V2",end.field = "V3")

IDRpeaks.H3K9ac.DS <- read.table("../01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_D_H3K9ac.IDR.uniq.narrowPeak",sep = "\t")%>%
  select(V1,V2,V3)%>%mutate(peak.id=str_c("DS_peak_",row.names(.)),width=V3-V2,treatment="DS")%>%unique()%>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,seqnames.field = "V1",start.field = "V2",end.field = "V3")

IDRpeaks.H3K9ac.RE <-read.table("../01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_RE_H3K9ac.IDR.uniq.narrowPeak",sep = "\t")%>%
  select(V1,V2,V3)%>%mutate(peak.id=str_c("RE_peak_",row.names(.)),width=V3-V2,treatment="RE")%>%unique()%>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,seqnames.field = "V1",start.field = "V2",end.field = "V3")

###Compare the relationships of the peaks under different conditions.
##Get the peak count of each class and plot it manually in AI.
#Get the total number of IDR peaks
length(IDRpeaks.H3K9ac.NC) # 18029
length(IDRpeaks.H3K9ac.DS) # 20000
length(IDRpeaks.H3K9ac.RE) # 19407

##I. DS vs NC relations:
#1. find total overlaps in all types and get the non-overlapping peak number
#NC peaks that can find overlaps with DS ones: 16241 (non-overlapping 1788)
findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.DS,maxgap=-1 ,type = "any")@from %>% unique()%>%length()
#DS peaks that can find overlaps with NC ones: 17381 (non-overlapping 2619)
findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.NC,maxgap=-1 ,type = "any")@from %>% unique()%>%length()

#2. find NC peaks with regions covering multiple DS peaks (DS peaks produced by split of NC peaks).
findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.NC,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #count in NC: 980

findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.NC,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums() #count in DS: 2211

#3. find DS peaks with regions covering multiple NC peaks.
findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.DS,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #Count in DS: 87

findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.DS,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums()  #Count in NC: 178

#4 The other overlapping peaks are classified into "other" group and not discussed in the MS.

#NC = 16241 - 980 - 178 = 15083
#DS = 17381 - 2211 - 87 = 15083

##II. RE vs DS relations:
#1. find total overlaps in all types and get the non-overlapping peak number
#DS peaks that can find overlaps with RE ones: 17807 (non-overlapping 2193)
findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.RE,maxgap=-1 ,type = "any")@from %>% unique()%>%length()
#RE peaks that can find overlaps with DS ones: 17401 (non-overlapping 2006)
findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.DS,maxgap=-1 ,type = "any")@from %>% unique()%>%length()

#2. find DS peaks with regions covering multiple RE peaks (RE peaks produced by split of DS peaks).
findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.DS,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #count in DS: 140

findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.DS,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums() #count in RE: 289

#3. find RE peaks with regions covering multiple DS peaks.
findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.RE,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #Count in RE: 484

findOverlaps(query = IDRpeaks.H3K9ac.DS ,subject = IDRpeaks.H3K9ac.RE,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums()  #Count in DS: 1039

#4 The other overlapping peaks are classified into "other" group and not discussed in the MS.

#DS = 17807 - 140 - 1039 = 16628
#RE = 17401 - 289 - 484 = 16628


##III. RE vs NC relations:
#1. find total overlaps in all types and get the non-overlapping peak number
#NC peaks that can find overlaps with RE ones: 17243 (non-overlapping 786)
findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.RE,maxgap=-1 ,type = "any")@from %>% unique()%>%length()
#RE peaks that can find overlaps with NC ones: 18053 (non-overlapping 1354)
findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.NC,maxgap=-1 ,type = "any")@from %>% unique()%>%length()

#2. find NC peaks with regions covering multiple RE peaks (RE peaks produced by split of NC peaks).
findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.NC,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #count in NC: 755

findOverlaps(query = IDRpeaks.H3K9ac.RE ,subject = IDRpeaks.H3K9ac.NC,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums() #count in RE: 1631

#3. find RE peaks with regions covering multiple NC peaks.
findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.RE,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%dim() #Count in RE: 64

findOverlaps(query = IDRpeaks.H3K9ac.NC ,subject = IDRpeaks.H3K9ac.RE,maxgap = -1, type = "any")@to %>%
  table%>%as.data.frame()%>%filter(Freq > 1)%>%select(Freq)%>%as.matrix%>%colSums()  #Count in DS: 130

#4 The other overlapping peaks are classified into "other" group and not discussed in the MS.

#NC = 17243 - 755 - 130 = 16358
#RE = 18053 - 1631 - 64 = 16358



####plot peak length (H3K9ac marked genomic span) under each condition
#Density plot
peak.width.all <- rbind(IDRpeaks.H3K9ac.NC%>%as.data.frame()%>%select(width,treatment),
                        IDRpeaks.H3K9ac.DS%>%as.data.frame()%>%select(width,treatment),
                        IDRpeaks.H3K9ac.RE%>%as.data.frame()%>%select(width,treatment))

peak.width.all$treatment <- factor(peak.width.all$treatment,levels = c("RE","DS","NC"),ordered = T)

pdf(file = "../02_plot/01_peak_span.distribute.pdf",height=2.28,width=3.72,bg=bg)
ggplot(peak.width.all,aes(x=log10(width),y=after_stat(count),color=treatment,group=treatment))+
  geom_density(size=0.75)  +
  scale_color_manual(values = condition.col)+
  #facet_wrap(condition~.,nrow = 1,ncol = 3,scales = "fixed")+
  xlab("Log10(H3K9ac marked region (bp))")+
  ylab("Peak count")
dev.off()

#Associated boxplot
pdf(file = "../02_plot/02_peak_span.boxplot.pdf",height=2.28,width=3.72,bg=bg)
ggplot(peak.width.all,aes(x=log10(width),y=treatment,fill=treatment))+
  geom_boxplot(outlier.alpha = 0.5,outlier.size = 0.3,width=0.8,outlier.fill = bg,outlier.colour = "grey50")  +
  scale_fill_manual(values = condition.col)+
  xlab("Log10(H3K9ac marked region (bp))")
dev.off()

#Statistic summary
aov(data = peak.width.all,formula = width~treatment)%>%TukeyHSD()

wilcox.test((peak.width.all%>%filter(treatment=="NC"))[,1]%>%log10,(peak.width.all%>%filter(treatment=="DS"))[,1]%>%log10) # < 2.2e-16
wilcox.test((peak.width.all%>%filter(treatment=="DS"))[,1]%>%log10,(peak.width.all%>%filter(treatment=="RE"))[,1]%>%log10) # < 2.2e-16
wilcox.test((peak.width.all%>%filter(treatment=="RE"))[,1]%>%log10,(peak.width.all%>%filter(treatment=="NC"))[,1]%>%log10) # < 2.2e-16


####=========================================================================================================######
####Plot the genomic coverage of the bined reads
####=========================================================================================================######

mh63.karyo <- read.table("~/genome_anno/MH63RS3/HM63RS3.rideogram.karyotype.txt",sep = "\t",header = T)

#The RIdeogram package only support the simultaneous plotting of two data sets in one command,
#So we plot the coverage by two steps and then combine the lines in AI.

#DS and NC
readCount.10kb%>%
  rowwise()%>%
  mutate(Value_1=max((MH_D_H3K9ac_1_IP.rmdup.bam+MH_D_H3K9ac_2_IP.rmdup.bam)/2-MH_input.rmdup.bam,0),
         Color_1="7B5141",

         Value_2=max((MH_N_H3K9ac_1_IP.rmdup.bam+MH_N_H3K9ac_2_IP.rmdup.bam)/2-MH_input.rmdup.bam,0),
         Color_2="86A47C",

         Value_3=max((MH_RE_H3K9ac_1_IP.rmdup.bam+MH_RE_H3K9ac_2_IP.rmdup.bam)/2-MH_input.rmdup.bam,0),
         Color_3="DFCEA1")%>%
  dplyr::rename(Chr=chr,Start=start,End=end)%>%
  select(Chr,Start,End,Value_1,Color_1,Value_2,Color_2)%>%
  as.data.frame()%>%
  RIdeogram::ideogram(
    karyotype = mh63.karyo,
    label = .,label_type = "line",
    output="../02_plot/ChIPseq_normalized_readsCoverage.ideo1.svg")

#RE
readCount.10kb%>%
  rowwise()%>%
  mutate(Value=max((MH_RE_H3K9ac_1_IP.rmdup.bam+MH_RE_H3K9ac_2_IP.rmdup.bam)/2-MH_input.rmdup.bam,0),
         Color="DFCEA1")%>%
  dplyr::rename(Chr=chr,Start=start,End=end)%>%
  select(Chr,Start,End,Value,Color)%>%
  as.data.frame()%>%
  RIdeogram::ideogram(
    karyotype = mh63.karyo,
    label = .,label_type = "line",
    output="../02_plot/ChIPseq_normalized_readsCoverage.ideo2.svg")



####=========================================================================================================######
####Annotate the genomic features associated with the IDR peaks
####=========================================================================================================######
anno.IDRpeaks.H3K9ac.NC <- IDRpeaks.H3K9ac.NC%>%
  annotatePeak(TxDb = mh63rs3.txdb,tssRegion = c(-1000,100),flankDistance=1000,addFlankGeneInfo=T)

anno.IDRpeaks.H3K9ac.DS <- IDRpeaks.H3K9ac.DS%>%
  annotatePeak(TxDb = mh63rs3.txdb,tssRegion = c(-1000,100),flankDistance=1000,addFlankGeneInfo=T)

anno.IDRpeaks.H3K9ac.RE <- IDRpeaks.H3K9ac.RE%>%
  annotatePeak(TxDb = mh63rs3.txdb,tssRegion = c(-1000,100),flankDistance=1000,addFlankGeneInfo=T)

#plot annoStat
pdf(file="../02_plot/01_IDRpeak.annoBar.pdf",height=2.28,width=3.72,bg=bg)
rbind(anno.IDRpeaks.H3K9ac.NC@annoStat%>%mutate(condition="NC"),
      anno.IDRpeaks.H3K9ac.DS@annoStat%>%mutate(condition="DS"),
      anno.IDRpeaks.H3K9ac.RE@annoStat%>%mutate(condition="RE"))%>%
  ggplot(data = ., aes(x=condition,fill=Feature,y=Frequency))+
  geom_bar(position='stack', stat='identity')+
  scale_fill_manual(values = ggsci::pal_uchicago()(9))
#scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 8,name = "Dark2"),ralcolors["RAL5012"]%>%as.vector()))
dev.off()

#Output annotated IDR peaks
anno.IDRpeaks.H3K9ac.NC@anno%>%as.data.frame%>%
  write.xlsx("../01_Peak_set/05_IDR_genomic_features/MH_N_H3K9ac.IDR.chipseeker.xlsx")

anno.IDRpeaks.H3K9ac.DS@anno%>%as.data.frame%>%
  write.xlsx("../01_Peak_set/05_IDR_genomic_features/MH_D_H3K9ac.IDR.chipseeker.xlsx")

anno.IDRpeaks.H3K9ac.RE@anno%>%as.data.frame%>%
  write.xlsx("../01_Peak_set/05_IDR_genomic_features/MH_RE_H3K9ac.IDR.chipseeker.xlsx")

####=========================================================================================================######
####Quantitative analysis on the differentially marked regions (DMR)
####=========================================================================================================######
#plot vocano plot to show DMR
voc.in <- rbind(anno_MH_D_H3K9ac_vs_MH_N_H3K9ac%>%as.data.frame()%>%select(-contains("Conc"))%>%mutate(cat="DS_vs_NC"),
                anno_MH_RE_H3K9ac_vs_MH_D_H3K9ac%>%as.data.frame()%>%select(-contains("Conc"))%>%mutate(cat="RE_vs_DS"),
                anno_MH_RE_H3K9ac_vs_MH_N_H3K9ac%>%as.data.frame()%>%select(-contains("Conc"))%>%mutate(cat="RE_vs_NC"))%>%
  select(Fold,FDR,annotation,geneId,cat)%>%
  merge(.,ricefungene,by="geneId",all.x=T)

voc.in <- voc.in%>%mutate(position=ifelse(str_detect(annotation,"Promoter")==T,"pro","other"))%>%
  mutate(dba=ifelse(Fold>1 & FDR < 0.05, "up",
                    ifelse(Fold < -1 & FDR < 0.05, "down","nc")))

#voc.lable <- c("SAPK4|OSPDK","SAPK10","SAPK1|OsSAPK1","RAB21|Rab16A|OsRab16A","OsSWI3C|OsCHB705|CHB705","OsSRO1c|BOC1","OsSIPP2C1|OsPP2C68|OsPP108","OsSAPK8|SAPK8","OsPP2C51",
#               "OsPP2C49","OsPP2C09|OsPP15|OsPYL2","OsNAC6|SNAC2","OsAREB8|OsAREB1|OsbZIP46|OsABF2|ABL1")
#voc.lable.data <- voc.in%>%filter(Gene_Symbol %in% voc.lable,dba != "nc",position=="pro")

pdf(file = "../02_plot/03_H3K9ac_DMR.volcano.pdf", width = 5.26,height =1.57,bg=bg)
ggplot() +

  geom_hex(data = voc.in%>%filter(dba == "nc") , aes(Fold, -log10(FDR)),bins=100)+
  scale_fill_gradient(low = "grey75" ,high = "black")+
  new_scale_fill()+

  geom_point(data = voc.in%>%filter(dba == "up",position=="other") , aes(Fold, -log10(FDR)),shape=23,size=0.9,color=bg,fill=alpha("orange",0.5))+
  geom_point(data = voc.in%>%filter(dba == "up",position=="pro") , aes(Fold, -log10(FDR)),shape=23,size=0.9,color=bg,fill=alpha("tomato",0.5))+

  geom_point(data = voc.in%>%filter(dba == "down",position=="other") , aes(Fold, -log10(FDR)),shape=23,size=0.9,color=bg,fill=alpha("#00bfff",0.5))+
  geom_point(data = voc.in%>%filter(dba == "down",position=="pro") , aes(Fold, -log10(FDR)),shape=23,size=0.9,color=bg,fill=alpha("#7fffd4",0.5))+
  #geom_point(data=voc.lable.data%>%filter(dba=="up") , aes(Fold, -log10(FDR)),shape=24,color="black",fill=alpha("#ff69b4",0.6),size=4)+
  #geom_point(data=voc.lable.data%>%filter(dba=="down") , aes(Fold, -log10(FDR)),shape=24,color="black",fill=alpha("#9370bd",0.6),size=4)+

  #scale_y_continuous(expand = c(0,0.2))+

  #geom_text_repel(data=voc.lable.data ,
  #                aes(Fold, -log10(FDR),label=Gene_Symbol),
  #                size=3,segment.size=0.2,show.legend = F,force=10)+

  facet_wrap(cat~.,scales = "free")

dev.off()

####=========================================================================================================######
####Map peak quantifications to genes
####=========================================================================================================######
#ChipSeeker package only annotate the peak to the nearest gene, other flanking genes are folded in the "flanking gene" column,
#therefore, many peak-overlaping genes beside the "nearest gene" will be lost in gene-based analysis if we use the peakAnno object directly.
#In this step we unfold the "flanking genes" whose promoter or body (flank_gene_distances within -1000 - 0) was covered by the peak,
#and map the peak intensity to them.
#Since non-promoter marked genes are hard to map and are not important in the following analysis, they are discarded in this step.
#DMG=differentially marked gene (by H3K9ac)
#MG=marked gene (by promoter H3K9ac)
#BA=binding aera
#DBA=differential binding aera
#-width,-strand,-annotation,geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds

BA.pro.D_vs_N.fullMG <- anno_MH_D_H3K9ac_vs_MH_N_H3K9ac%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

BA.pro.RE_vs_D.fullMG  <- anno_MH_RE_H3K9ac_vs_MH_D_H3K9ac%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

BA.pro.RE_vs_N.fullMG  <- anno_MH_RE_H3K9ac_vs_MH_N_H3K9ac%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

#Output the results:
write.xlsx(file = "../08_Promoter_marked_genes/DBA.DS_VS_NC.AllPromoterMarkedGenes.xlsx",x=BA.pro.D_vs_N.fullMG)
write.xlsx(file = "../08_Promoter_marked_genes/DBA.RE_VS_DS.AllPromoterMarkedGenes.xlsx",x=BA.pro.RE_vs_D.fullMG)
write.xlsx(file = "../08_Promoter_marked_genes/DBA.RE_VS_NC.AllPromoterMarkedGenes.xlsx",x=BA.pro.RE_vs_N.fullMG)

####=========================================================================================================######
####Correlation analysis with expression data.
####=========================================================================================================######

#Use the log2FC of BAs to correlate the expression log2FC of MGs in the DBAs.

#DS vd NC
deg.dn.FC <- read.delim("../../../03_RNAseq/02_DEG/01_DEG_table/allGenes_MH_D_vs_MH_N.txt",sep = "\t",header = T)%>%
  select(Row.names,log2FoldChange)%>%rename(Expression_Fold_Change=log2FoldChange)

h3k9ac.dn.FC <- DBA.pro.D_vs_N.fullMG%>%select(geneId,Fold)%>%rename(H3K9ac_Fold_Change=Fold)

h3k9ac.exp.cor.in <- merge(deg.dn.FC,h3k9ac.dn.FC,by.x="Row.names",by.y="geneId",all=F)%>%na.omit()

pdf("../02_plot/04_correlation.H3K9ac_Exp.DN.pdf",height = 1.75,width = 2.75,bg=bg)
ggplot(data = h3k9ac.exp.cor.in,aes(x=Expression_Fold_Change,y = H3K9ac_Fold_Change))+
  stat_bin_hex(bins = 100,size=5)+
  scale_fill_gradientn(colours = viridis::plasma(255))+
  #scale_fill_gradient2(high = "tomato",low = "#7b68ee",mid = "#40e0d0",midpoint = 125)+
  xlab("Expression Fold Change(DS vs NC)")+
  ylab("H3K9ac Fold Change(DS vs NC)")+
  geom_smooth(method = "lm",se = T,color=alpha("grey25",0.4))
dev.off()

cor.test(h3k9ac.exp.cor.in$Expression_Fold_Change,h3k9ac.exp.cor.in$H3K9ac_Fold_Change,method = "pearson") #PPC=0.6314239, p-value < 2.2e-16

#RE vs DS
deg.rd.FC <- read.delim("../../../03_RNAseq/02_DEG/01_DEG_table/allGenes_MH_RE_vs_MH_D.txt",sep = "\t",header = T)%>%
  select(Row.names,log2FoldChange)%>%rename(Expression_Fold_Change=log2FoldChange)

h3k9ac.rd.FC <- DBA.pro.RE_vs_D.fullMG%>%select(geneId,Fold)%>%rename(H3K9ac_Fold_Change=Fold)

h3k9ac.exp.cor.rd.in <- merge(deg.rd.FC,h3k9ac.rd.FC,by.x="Row.names",by.y="geneId",all=F)%>%na.omit()

pdf("../02_plot/05_correlation.H3K9ac_Exp.RD.pdf",height = 1.75,width = 2.75,bg=bg)
ggplot(data = h3k9ac.exp.cor.rd.in,aes(x=Expression_Fold_Change,y = H3K9ac_Fold_Change))+
  stat_bin_hex(bins = 100,size=5)+
  scale_fill_gradientn(colours = viridis::plasma(255))+
  #scale_fill_gradient2(high = "tomato",low = "#7b68ee",mid = "#40e0d0",midpoint = 200)+
  xlab("Expression Fold Change(RE vs DS)")+
  ylab("H3K9ac Fold Change(RE vs DS)")+
  geom_smooth(method = "lm",se = T,color=alpha("grey25",0.4))
dev.off()

cor.test(h3k9ac.exp.cor.rd.in$Expression_Fold_Change,h3k9ac.exp.cor.rd.in$H3K9ac_Fold_Change,method = "pearson") #PPC=0.6634867, p-value < 2.2e-16

#RE vs NC
deg.rn.FC <- read.delim("../../../03_RNAseq/02_DEG/01_DEG_table/allGenes_MH_RE_vs_MH_N.txt",sep = "\t",header = T)%>%
  select(Row.names,log2FoldChange)%>%rename(Expression_Fold_Change=log2FoldChange)

h3k9ac.rn.FC <- DBA.pro.RE_vs_N.fullMG%>%select(geneId,Fold)%>%rename(H3K9ac_Fold_Change=Fold)

h3k9ac.exp.cor.rn.in <- merge(deg.rn.FC,h3k9ac.rn.FC,by.x="Row.names",by.y="geneId",all=F)%>%na.omit()

pdf("../02_plot/05_correlation.H3K9ac_Exp.RN.pdf",height = 1.75,width = 2.75,bg=bg)
ggplot(data = h3k9ac.exp.cor.rn.in,aes(x=Expression_Fold_Change,y = H3K9ac_Fold_Change))+
  stat_bin_hex(bins = 100,size=5)+
  scale_fill_gradientn(colours = viridis::plasma(255))+
  #scale_fill_gradient2(high = "tomato",low = "#7b68ee",mid = "#40e0d0",midpoint = 200)+
  xlab("Expression Fold Change(RE vs NC)")+
  ylab("H3K9ac Fold Change(RE vs NC)")+
  geom_smooth(method = "lm",se = T,color=alpha("grey25",0.4))
dev.off()

cor.test(h3k9ac.exp.cor.rn.in$Expression_Fold_Change,h3k9ac.exp.cor.rn.in$H3K9ac_Fold_Change,method = "pearson") #PPC=0.1578765, p-value < 2.2e-16



####=========================================================================================================######
####Get the full list of promoter-marked genes by IDR peaks
####=========================================================================================================######

#Similar to the full list of BA-associated genes, we get all genes with promoter associated/covered by IDR peak.

IDRpeaks.pro.DS.fullMG <- anno.IDRpeaks.H3K9ac.DS%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

IDRpeaks.pro.NC.fullMG <- anno.IDRpeaks.H3K9ac.NC%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

IDRpeaks.pro.RE.fullMG <- anno.IDRpeaks.H3K9ac.RE%>%as.data.frame()%>%filter(annotation=="Promoter")%>%
  select(-width,-strand,-annotation,-geneId,-geneId,-transcriptId,-distanceToTSS,-geneChr,-geneStart,-geneEnd,-geneLength,-geneStrand,-flank_txIds)%>%
  tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep = ";",convert = T)%>%as.data.frame()%>%unique()%>%
  filter(flank_gene_distances > -1000, flank_gene_distances <=0)%>%rename(geneId=flank_geneIds)%>%
  left_join(.,ricefungene,by="geneId")

#Output the results:
write.xlsx(file = "../08_Promoter_marked_genes/IDRpeaks.NC.AllPromoterMarkedGenes.xlsx",x=IDRpeaks.pro.NC.fullMG)
write.xlsx(file = "../08_Promoter_marked_genes/IDRpeaks.DS.AllPromoterMarkedGenes.xlsx",x=IDRpeaks.pro.DS.fullMG)
write.xlsx(file = "../08_Promoter_marked_genes/IDRpeaks.RE.AllPromoterMarkedGenes.xlsx",x=IDRpeaks.pro.RE.fullMG)
