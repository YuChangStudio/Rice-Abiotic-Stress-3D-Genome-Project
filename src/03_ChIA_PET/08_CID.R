####=========================================================================================================######
##CID: genomic regions referring to continuous local loops connected by shared anchors
##loops with shared anchors may be regulated by identical trans-factors binding to\flanking the anchors,
##thus genes connected by the loops are more likely to be co-regulated.
####=========================================================================================================######

####=========================================================================================================######
####Load overall configs
####=========================================================================================================######
library(multcomp)
library(ggplot2)
library(GenomicRanges)
library(ral)
library(stringr)
library(dplyr)

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
      #legend.position=c(0.056,0.88),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white",colour = "black"),
      legend.key = element_blank(),
      plot.background = element_blank())%>%theme_set()


#color
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"
col.bzip23 <- ralcolors["RAL6027"]%>%as.vector #"#81C0BB"

condition.col<- c("NC"=col.NC,"DS"=col.DS,"RE"=col.RE,"bzip23"=col.bzip23)

#TPM (Tag per Million in this case)
get_TPM <- function(counts, lengths) {
  stopifnot(is.numeric(counts))
  stopifnot(is.numeric(lengths))
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

#Data header
loacl.loops.header <- c("chrom1","start1",	"end1",	"chrom2",	"start2",	"end2","loop.index","ipet.counts")

####=========================================================================================================######
####load local loops
####=========================================================================================================######
##load loacl loops (loop span <= 200,000 bp) identified in the loop profile analysis.
#the loop index is in accordance with the full loop set.
DS.local.loops <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/D.local.loop.bedpe",header = F)
NC.local.loops <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/N.local.loop.bedpe",header = F)
RE.local.loops <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/RE.local.loop.bedpe",header = F)
bzip23.local.loops <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/bzip23.local.loop.bedpe",header = F)


colnames(NC.local.loops) <- loacl.loops.header
colnames(DS.local.loops) <- loacl.loops.header
colnames(RE.local.loops) <- loacl.loops.header
colnames(bzip23.local.loops) <- loacl.loops.header

####=========================================================================================================######
####Find continuously connecting loops
####=========================================================================================================######
#Use bedtools pairtopair to find loops with at least one shared anchor.
#Extract those regions (pairwise loop linkage, PLL) to bedgraph format.
library(bedtoolsr)
options(bedtools.path = "/Users/juliuschiong/opt/anaconda3/bin/")

for (x in c("DS","NC","RE","bzip23")){

  loops = get(paste0(x,".local.loops"))

  bt.pairtopair(a=loops,b=loops,is = T,rdn = T,type = "either")%>%
    merge(.,loops%>%select(loop.index,ipet.counts),by.x="V7",by.y="loop.index",all.x=T)%>%
    merge(.,loops%>%select(loop.index,ipet.counts),by.x="V15",by.y="loop.index",all.x=T)%>%
    rowwise%>%
    mutate(loop1.index=str_sort(c(V7,V15),numeric = T)[1],loop2.index=str_sort(c(V7,V15),numeric = T)[2])%>%
    mutate(chr=V1,
           str=unlist(min(V2,V3,V5,V6,V10,V11,V13,V14)),
           end=unlist(max(V2,V3,V5,V6,V10,V11,V13,V14)),
           ipet=ipet.counts.x+ipet.counts.y)%>%
    select(chr,str,end,ipet,loop1.index,loop2.index)%>%
    as.data.frame()%>%
    arrange(chr,str,end,loop1.index,loop2.index)%>%
    unique()%>%
    mutate(width=end-str)%>%
    assign(paste0(x,".PLL"),.,envir = .GlobalEnv)

  rm(loops)
}

rm(x,loops)

####=========================================================================================================######
####Merge the PLL regions to form candidate loop clusters (preDomains)
####=========================================================================================================######

#merge PLL regions by bedtools merge

for (x in c("DS","NC","RE","bzip23")) {
  PLL=get(paste0(x,".PLL"))
  bt.merge(i = PLL,c = paste0(4,",",5,",",6),o = paste0("sum",",","distinct",",","distinct") ,delim = "\\;")%>%
    rowwise()%>%
    mutate(loop.index=unique(c(str_split(V5,pattern = ";")%>%unlist(),
                               str_split(V6,pattern = ";")%>%unlist()))%>%
             paste(collapse = ";"))%>%
    mutate(uniq.loop.count=length(str_split(loop.index,pattern = ";")%>%unlist()))%>%
    mutate(width=V3-V2)%>%
    dplyr::rename(chr=V1,start=V2,end=V3,total.ipet=V4)%>%
    select(chr,start,end,total.ipet,uniq.loop.count,width,loop.index)%>%
    as.data.frame()%>%
    assign(paste0(x,".preDomain"),.,envir = .GlobalEnv)
}

rm(x,PLL)

####=========================================================================================================######
####Filter insignificant pre-domains by loop count
####=========================================================================================================######

quantile(DS.preDomain$total.ipet,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(DS.preDomain$uniq.loop.count,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(DS.preDomain$width,probs=c(0.1,0.25,0.50,0.75,0.90,1))

quantile(NC.preDomain$total.ipet,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(NC.preDomain$uniq.loop.count,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(NC.preDomain$width,probs=c(0.1,0.25,0.50,0.75,0.90,1))

quantile(RE.preDomain$total.ipet,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(RE.preDomain$uniq.loop.count,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(RE.preDomain$width,probs=c(0.1,0.25,0.50,0.75,0.90,1))

quantile(bzip23.preDomain$total.ipet,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(bzip23.preDomain$uniq.loop.count,probs=c(0.1,0.25,0.50,0.75,0.90,1))
quantile(bzip23.preDomain$width,probs=c(0.1,0.25,0.50,0.75,0.90,1))

#Discard the pre-domains with loop count=2 (the 10% quantile)
#The rest of the pre-domain regions are regarded as Chromatin Interacting Domains (CIDs)

DS.CID <- DS.preDomain%>%filter(uniq.loop.count > 2)
NC.CID <- NC.preDomain%>%filter(uniq.loop.count > 2)
RE.CID <- RE.preDomain%>%filter(uniq.loop.count > 2)
bzip23.CID <- bzip23.preDomain%>%filter(uniq.loop.count > 2)

#output the regions

for (x in c("DS","NC","RE","bzip23")) {
  CID=get(paste0(x,".CID"))

  #output as standard bedgraph-6 format
  #use uniq.loop.count as score

  CID%>%
    mutate(name=paste0(x,"_CID_",rownames(.)),score=uniq.loop.count,strand="*")%>%
    select(chr,start,end,name,score,strand)%>%
    write.table(file = paste0("../01_CID_regions/",x,"_CID.bed"),.,sep = "\t",col.names = F,row.names = F,quote = F)

  CID%>%
    write.table(file = paste0("../01_CID_regions/",x,"_CID.fullInfo.txt"),.,sep = "\t",col.names = T,row.names = F,quote = F)
}

rm(x,CID)

####==================================Analysis on the CID regions============================================######
####=========================================================================================================######
####Compare the CID covered regions
####=========================================================================================================######

##Anova analysis
CID.span <- rbind(NC.CID%>%mutate(condition="NC",span=end-start)%>%select(condition,span),
                  DS.CID%>%mutate(condition="DS",span=end-start)%>%select(condition,span),
                  RE.CID%>%mutate(condition="RE",span=end-start)%>%select(condition,span),
                  bzip23.CID%>%mutate(condition="bzip23",span=end-start)%>%select(condition,span))%>%
  mutate(condition=factor(condition,levels = c("NC","DS","RE","bzip23"),ordered = T))

CID.span%>%mutate(log2.span=log2(span))%>%lm(log2.span~condition,data = .)%>%aov()%>%TukeyHSD()

CID.span.hsd <- CID.span%>%mutate(log2.span=log2(span))%>%lm(log2.span~condition,data = .)%>%
  glht(.,linfct = mcp(condition = "Tukey"))%>%cld(decreasing=T)

CID.span.hsd<-CID.span.hsd$mcletters$Letters%>%as.data.frame()%>%rename(hsd=".")%>%mutate(condition=rownames(.))

##plot
##For MH63 comparision

CID.span.MH63 <- CID.span%>%filter(condition != "bzip23")%>%mutate(condition=factor(condition,levels = c("NC","DS","RE"),ordered = T))

pdf("../02_CID_profile_analysis/CID.spans.MH63.pdf",height = 3,width = 3.2,bg=bg)
ggplot(data =CID.span.MH63 ,aes(x=condition,y=log2(span)))+
  geom_violin(aes(fill=condition),trim = F,scale = "count")+
  geom_boxplot(fill="white",width=0.3,outlier.shape = NA)+
  stat_summary(fun.y = mean,size=0.1)+
  scale_fill_manual(values = condition.col)+
  geom_text(data=CID.span.hsd[1:3,],aes(x=condition,y=1.1*log2(max(CID.span.MH63$span)),label=hsd),inherit.aes = F)+
  #ggsignif::geom_signif(comparisons = list(c("DS","NC"),c("RE","NC"),c("RE","DS")),step_increase = 0.15)+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab("log2(CID span)")
dev.off()

##For bizp23 comparision
pdf("../02_CID_profile_analysis/CID.spans.MH63_and_bzip23.pdf",height = 3,width = 3.5,bg=bg)
ggplot(data = CID.span ,aes(x=condition,y=log2(span)))+
  geom_violin(aes(fill=condition))+
  geom_boxplot(fill="white",width=0.3,outlier.shape = NA)+
  stat_summary(fun.y = mean,size=0.1)+
  scale_fill_manual(values = condition.col)+
  geom_text(data=CID.span.hsd,aes(x=condition,y=1.05*log2(max(CID.span$span)),label=hsd),inherit.aes = F)+
  #ggsignif::geom_signif(comparisons = list(c("DS","NC"),c("RE","NC"),c("RE","DS")),step_increase = 0.15)+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab("log2(CID span)")
dev.off()

#Get the total CID coverage under each condtion
CID.span%>%group_by(condition)%>%mutate(sum=sum(span))%>%select(condition,sum)%>%unique()
CID.span%>%group_by(condition)%>%mutate(sum=n())%>%select(condition,sum)%>%unique()
CID.span%>%group_by(condition)%>%mutate(sum=mean(span))%>%select(condition,sum)%>%unique()
####=========================================================================================================######
####Compare gene expressions in/out the CIDs
####=========================================================================================================######
library(GenomicRanges)
###Get gene expression levels in the CIDs (all genes are included).
mh63.gene.bed <- read.table("~/genome_anno/MH63RS3/MH63RS3.gene.bed")[,1:6]%>%unique()

TPM.1 <- read.table(file = "../../../03_RNAseq_MH63/01_data_table_and_sample_list/gene.tpm",header = T)%>%
  mutate(DS=(MH_D_1+MH_D_2)/2,NC=(MH_N_1+MH_N_2)/2,RE=(MH_RE_1+MH_RE_2)/2)%>%
  select(gene_id,DS,NC,RE)%>%reshape::melt()%>%dplyr::rename(condition=variable,TPM=value)

TPM.2 <- read.table(file = "../../../06_RNAseq_bzip23/00_quant_matrix/gene.tpm",header = T)%>%
  mutate(bzip23=(bZIP23_D1+bZIP23_D2)/2)%>%
  select(gene_id,bzip23)%>%reshape::melt()%>%dplyr::rename(condition=variable,TPM=value)

TPM.all <- rbind(TPM.1,TPM.2)

#Genes in DS CIDs
DS.CID.gene <- findOverlaps(query  = mh63.gene.bed%>%makeGRangesFromDataFrame(seqnames.field  = "V1",start.field = "V2",end.field = "V3"),
                            subject  = DS.CID%>%makeGRangesFromDataFrame(seqnames.field  = "chr",start.field = "start",end.field = "end")
                            ,type = "within")@from%>%unique()

DS.CID.gene <- mh63.gene.bed[DS.CID.gene,]$V4

length(DS.CID.gene) #14076

#Genes in NC CIDs
NC.CID.gene <- findOverlaps(query  = mh63.gene.bed%>%makeGRangesFromDataFrame(seqnames.field  = "V1",start.field = "V2",end.field = "V3"),
                            subject  = NC.CID%>%makeGRangesFromDataFrame(seqnames.field  = "chr",start.field = "start",end.field = "end")
                            ,type = "within")@from%>%unique()

NC.CID.gene <- mh63.gene.bed[NC.CID.gene,]$V4

length(NC.CID.gene) #23386

#Genes in RE CIDs
RE.CID.gene <- findOverlaps(query  = mh63.gene.bed%>%makeGRangesFromDataFrame(seqnames.field  = "V1",start.field = "V2",end.field = "V3"),
                            subject  = RE.CID%>%makeGRangesFromDataFrame(seqnames.field  = "chr",start.field = "start",end.field = "end")
                            ,type = "within")@from%>%unique()

RE.CID.gene <- mh63.gene.bed[RE.CID.gene,]$V4

length(RE.CID.gene) #12478

#Genes in bzip23 CIDs
bzip23.CID.gene <- findOverlaps(query  = mh63.gene.bed%>%makeGRangesFromDataFrame(seqnames.field  = "V1",start.field = "V2",end.field = "V3"),
                                subject  = bzip23.CID%>%makeGRangesFromDataFrame(seqnames.field  = "chr",start.field = "start",end.field = "end")
                                ,type = "within")@from%>%unique()

bzip23.CID.gene <- mh63.gene.bed[bzip23.CID.gene,]$V4

length(bzip23.CID.gene) #15424


#compare the expressions of the three categories of genes
pdf("../02_CID_profile_analysis/CID.geneExp.pdf",height = 3,width = 9,bg=bg)
rbind(TPM.all%>%filter(gene_id %in% DS.CID.gene,condition== "DS")%>%mutate(cat="in_CIDs"),
      TPM.all%>%filter(!(gene_id %in% DS.CID.gene),condition== "DS")%>%mutate(cat="out_CIDs"),
      TPM.all%>%filter(gene_id %in% NC.CID.gene,condition== "NC")%>%mutate(cat="in_CIDs"),
      TPM.all%>%filter(!(gene_id %in% NC.CID.gene),condition== "NC")%>%mutate(cat="out_CIDs"),
      TPM.all%>%filter(gene_id %in% RE.CID.gene,condition== "RE")%>%mutate(cat="in_CIDs"),
      TPM.all%>%filter(!(gene_id %in% RE.CID.gene),condition== "RE")%>%mutate(cat="out_CIDs"),
      TPM.all%>%filter(gene_id %in% bzip23.CID.gene,condition== "bzip23")%>%mutate(cat="in_CIDs"),
      TPM.all%>%filter(!(gene_id %in% bzip23.CID.gene),condition== "bzip23")%>%mutate(cat="out_CIDs"))%>%
  filter(TPM>1)%>%
  ggplot(data = .,aes(x=cat,y=log2(TPM)))+
  geom_violin(aes(fill=condition),trim = F,scale = "count")+
  geom_boxplot(fill="white",width=0.3,outlier.shape = NA)+
  stat_summary(fun.y = mean,size=0.1)+
  scale_fill_manual(values = condition.col)+
  ggsignif::geom_signif(comparisons = list(c("in_CIDs","out_CIDs")),step_increase = 0.15)+
  facet_wrap(.~condition,ncol = 4)+
  #theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab("log2(TPM of expressed genes)")
dev.off()

####=========================================================================================================######
####Compare the position relationships of CIDs under each condition
####=========================================================================================================######

#=================DS vs NC
#get CID counts
NC.CID.grange <- makeGRangesFromDataFrame(NC.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
DS.CID.grange <- makeGRangesFromDataFrame(DS.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
length(NC.CID.grange) #766
length(DS.CID.grange) #640

#find total overlaps in all types
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "any")@from%>%unique()%>%length() #count in NC:476
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "any")@to%>%unique()%>%length() #count in DS:624

#find equal CIDs.
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "equal")@from%>%unique()%>%length() #count in NC:38
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "equal")@to%>%unique()%>%length() #count in DS:38

#find NC CIDs that are inside DS CID spans
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in NC:78-38=40
findOverlaps(query  = NC.CID.grange,subject  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in DS:76-38=38

#find DS CIDs that are inside NC CID spans
findOverlaps(subject  = NC.CID.grange,query  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in DS:503-38=465
findOverlaps(subject  = NC.CID.grange,query  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in NC:362-38=324

#None-overlap equals total-type.any
#None-overlap.NC=766-476=290
#None-overlap.DS=640-624=16

#other types of interactions equals type.any-type.within-type.harbor-equal
#other.NC=476-40-324-38=74
#other.DS=624-38-465-38=83



#=================RE vs DS
#get CID counts
RE.CID.grange <- makeGRangesFromDataFrame(RE.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
DS.CID.grange <- makeGRangesFromDataFrame(DS.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
length(RE.CID.grange) #750
length(DS.CID.grange) #640

#find total overlaps in all types
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "any")@from%>%unique()%>%length() #count in RE:625
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "any")@to%>%unique()%>%length() #count in DS:532

#find equal CIDs.
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "equal")@from%>%unique()%>%length() #count in RE:88
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "equal")@to%>%unique()%>%length() #count in DS:88

#find RE CIDs that are inside DS CID spans
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in RE:401-88=313
findOverlaps(query  = RE.CID.grange,subject  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in DS:331-88=243

#find DS CIDs that are inside RE CID spans
findOverlaps(subject  = RE.CID.grange,query  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in DS:221-88=133
findOverlaps(subject  = RE.CID.grange,query  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in RE:209-88=121

#None-overlap equals total-type.any
#None-overlap.RE=750-625=125
#None-overlap.DS=640-532=108

#other types of interactions equals type.any-type.within-type.harbor-equal
#other.RE=625-313-121-88=103
#other.DS=532-243-133-88=68

#=================RE vs NC
#get CID counts
RE.CID.grange <- makeGRangesFromDataFrame(RE.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
NC.CID.grange <- makeGRangesFromDataFrame(NC.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
length(RE.CID.grange) #750
length(NC.CID.grange) #766

#find total overlaps in all types
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "any")@from%>%unique()%>%length() #count in RE:741
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "any")@to%>%unique()%>%length() #count in NC:487

#find equal CIDs.
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "equal")@from%>%unique()%>%length() #count in RE:55
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "equal")@to%>%unique()%>%length() #count in NC:55

#find RE CIDs that are inside NC CID spans
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "within")@from%>%unique()%>%length() #count in RE:692-55=637
findOverlaps(query  = RE.CID.grange,subject  = NC.CID.grange,type = "within")@to%>%unique()%>%length() #count in NC:453-55=398

#find NC CIDs that are inside RE CID spans
findOverlaps(subject  = RE.CID.grange,query  = NC.CID.grange,type = "within")@from%>%unique()%>%length() #count in NC:75-55=20
findOverlaps(subject  = RE.CID.grange,query  = NC.CID.grange,type = "within")@to%>%unique()%>%length() #count in RE:75-55=20

#None-overlap equals total-type.any
#None-overlap.RE=750-741=9
#None-overlap.NC=766-487=279

#other types of interactions equals type.any-type.within-type.harbor-equal
#other.RE=741-637-20-55=29
#other.NC=487-398-20-55=14

#=================DS vs bzip23
#get CID counts
bzip23.CID.grange <- makeGRangesFromDataFrame(bzip23.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
DS.CID.grange <- makeGRangesFromDataFrame(DS.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
length(bzip23.CID.grange) #439
length(DS.CID.grange) #640

#find total overlaps in all types
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "any")@from%>%unique()%>%length() #count in bzip23:405
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "any")@to%>%unique()%>%length() #count in DS:484

#find equal CIDs.
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "equal")@from%>%unique()%>%length() #count in bzip23:43
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "equal")@to%>%unique()%>%length() #count in DS:43

#find bzip23 CIDs that are inside DS CID spans
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in bzip23:154-43=111
findOverlaps(query  = bzip23.CID.grange,subject  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in DS:149-43=106

#find DS CIDs that are inside bzip23 CID spans
findOverlaps(subject  = bzip23.CID.grange,query  = DS.CID.grange,type = "within")@from%>%unique()%>%length() #count in DS:263-43=220
findOverlaps(subject  = bzip23.CID.grange,query  = DS.CID.grange,type = "within")@to%>%unique()%>%length() #count in bzip23:209-43=166

#None-overlap equals total-type.any
#None-overlap.bzip23=439-405=34
#None-overlap.DS=640-484=156

#other types of interactions equals type.any-type.within-type.harbor-equal
#other.bzip23=405-111-166-43=85
#other.DS=484-106-220-43=115

#=================NC vs bzip23
#get CID counts
bzip23.CID.grange <- makeGRangesFromDataFrame(bzip23.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
NC.CID.grange <- makeGRangesFromDataFrame(NC.CID,seqnames.field  = "chr",start.field = "start",end.field = "end")%>%unique()
length(bzip23.CID.grange) #439
length(NC.CID.grange) #766

#find total overlaps in all types
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "any")@from%>%unique()%>%length() #count in bzip23:435
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "any")@to%>%unique()%>%length() #count in NC:384

#find equal CIDs.
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "equal")@from%>%unique()%>%length() #count in bzip23:21
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "equal")@to%>%unique()%>%length() #count in NC:21

#find bzip23 CIDs that are inside NC CID spans
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "within")@from%>%unique()%>%length() #count in bzip23:313-21=292
findOverlaps(query  = bzip23.CID.grange,subject  = NC.CID.grange,type = "within")@to%>%unique()%>%length() #count in NC:249-21=228

#find NC CIDs that are inside bzip23 CID spans
findOverlaps(subject  = bzip23.CID.grange,query  = NC.CID.grange,type = "within")@from%>%unique()%>%length() #count in NC:76-21=55
findOverlaps(subject  = bzip23.CID.grange,query  = NC.CID.grange,type = "within")@to%>%unique()%>%length() #count in bzip23:69-21=48

#None-overlap equals total-type.any
#None-overlap.bzip23=439-435=4
#None-overlap.NC=766-384=382

#other types of interactions equals type.any-type.within-type.harbor-equal
#other.bzip23=435-292-48-21=74
#other.NC=384-228-55-21=80

#================================================
###The jaccard values will be calculated using bedtools jaccard in shell.
###Relationships were plotted by AI (CID_relations.ai).

####=========================================================================================================######
####Plot CID ideogram
####=========================================================================================================######

mh63.karyo <- read.table("~/genome_anno/MH63RS3/HM63RS3.rideogram.karyotype.txt",sep = "\t",header = T)

NC.CID%>%mutate(V4=1)%>%dplyr::rename(Chr=chr,Start=start,End=end,Value=V4)%>%
  RIdeogram::ideogram(karyotype = mh63.karyo,overlaid=.,colorset1=c(col.NC),output="../02_CID_profile_analysis/NC.CID.ideo.svg")

DS.CID%>%mutate(V4=1)%>%dplyr::rename(Chr=chr,Start=start,End=end,Value=V4)%>%
  RIdeogram::ideogram(karyotype = mh63.karyo,overlaid=.,colorset1=c(col.DS),output="../02_CID_profile_analysis/DS.CID.ideo.svg")

RE.CID%>%mutate(V4=1)%>%dplyr::rename(Chr=chr,Start=start,End=end,Value=V4)%>%
  RIdeogram::ideogram(karyotype = mh63.karyo,overlaid=.,colorset1=c(col.RE),output="../02_CID_profile_analysis/RE.CID.ideo.svg")

bzip23.CID%>%mutate(V4=1)%>%dplyr::rename(Chr=chr,Start=start,End=end,Value=V4)%>%
  RIdeogram::ideogram(karyotype = mh63.karyo,overlaid=.,colorset1=c(col.bzip23),output="../02_CID_profile_analysis/bzip23.CID.ideo.svg")
