####=========================================================================================================######
####Load overall configs
####=========================================================================================================######

library(ChIPseeker)
library(ggplot2)
library(openxlsx)
library(GenomicRanges)
library(parallel)
library(openxlsx)
library(stringr)
library(ggsignif)
library(ggplot2)
library(ral)
library(ggsignif)
library(multcomp)
library(dplyr)

#color
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"
col.bzip23 <- ralcolors["RAL6027"]%>%as.vector #"#81C0BB"

condition.col<- c("NC"=col.NC,"DS"=col.DS,"RE"=col.RE,"bzip23"=col.bzip23)

#ggplot2 theme
theme(panel.background = element_rect(fill = "transparent",colour='black',size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0.4,"lines"),
      axis.ticks = element_line(color='black',size =0.45 ),
      axis.line = element_line(colour = "black",size = 0),
      axis.title.x=element_text(colour='black', size=12,face = "bold"),
      axis.title.y=element_text(colour='black', size=12,face = "bold",angle = 90),
      axis.text.x=element_text(colour='black',size=8),
      axis.text.y = element_text(color = "black",size = 8),
      strip.text.x = element_text(angle = 0,size = 8,face = "bold"),
      #legend.position=c(0.056,0.88),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white",colour = "black"),
      legend.key = element_blank(),
      plot.background = element_blank())%>%theme_set()


##Functional gene annotations derived from funRiceGenes (March, 2021).
mh63.encode <- read.xlsx("~/genome_anno/MH63RS3/MH63RS3_encode_Sep2022.xlsx")

##Chr length
mh63.chrlength <- read.table("~/genome_anno/MH63RS3/HM63RS3.rideogram.karyotype.txt",sep = "\t",header = T)%>%
  select(Chr,End)

##Expression data
#tpm1: MH63 series:
tpm.1 <- read.table("../../03_RNAseq_MH63/01_data_table_and_sample_list/gene.tpm",header = T,sep = "\t")%>%
  mutate(DS=(MH_D_1+MH_D_2)/2,NC=(MH_N_1+MH_N_2)/2,RE=(MH_RE_1+MH_RE_2)/2)%>%
  #filter(expression_under_DS>1 | expression_under_NC>1 | expression_under_RE >1)%>%
  select(gene_id,DS,NC,RE)
#reshape2::melt()%>%
#rename(TPM=value,condition=variable)

#tpm2: bzip23 series:
tpm.2 <- read.table("../../06_RNAseq_bzip23/00_quant_matrix/gene.tpm",header = T,sep = "\t")%>%
  mutate(bzip23=(bZIP23_D1+bZIP23_D2)/2)%>%
  select(gene_id,bzip23)
#reshape2::melt()%>%
#rename(TPM=value,condition=variable)

tpm.all <- left_join(tpm.1,tpm.2)%>% #Get full quantification list.
  rename(TPM.DS=DS,TPM.NC=NC,TPM.RE=RE,TPM.bzip23=bzip23)

rm(tpm.1,tpm.2)

####=========================================================================================================######
#Filter iPET clusters and optimize the iPET count and P-value cut-off
####=========================================================================================================######
sample.tag <- c("N","D","RE","bzip23")

header.cluster <- c("chrom1","start1",	"end1",	"chrom2",	"start2",	"end2",
                    "ipet.counts",	"type",	"distance",	"tag.count.within.anchor.1",
                    "tag.count.within.anchor.2",	"p-value",	"p.adjust",	"-log10(p-value)",	"-log10(p.adjust)")

lapply(sample.tag, function(x){

  prefix <- "../../02_ChIA-PET/01_ChIA-PET_tool_out/H3K9ac_"
  file <- paste0(prefix,x,"/03_clusters/H3K9ac_",x,"_pooled.cluster.FDRfiltered.txt")

  df <- read.delim(paste0(file),header=F)

  colnames(df) <- header.cluster

  assign(paste0("cluster.",x),df,envir = .GlobalEnv)

})


#Filter confident PET clusters (or so called loops) at ipet.count >=6, FDR<0.05, and distance > 0

cluster.N%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==1,distance > 0)%>%dim() #19439
cluster.N%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==0)%>%dim() #4796

cluster.D%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==1,distance > 0)%>%dim() #9262
cluster.D%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==0)%>%dim() #315

cluster.RE%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==1,distance > 0)%>%dim() #8687
cluster.RE%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==0)%>%dim() #1939

cluster.bzip23%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==1,distance > 0)%>%dim() #8325
cluster.bzip23%>%filter(ipet.counts >= 6 , p.adjust < 0.05,type==0)%>%dim() #41

####=========================================================================================================######
#Load loop (iPET clusters with iPET count >= 6) results based on the ChIP-Seq peaks.
####=========================================================================================================######

#Cluster.txt header
header.loop <- c("chrom1","start1",	"end1",	"chrom2",	"start2",	"end2","ipet.counts","distance","p.adjust","type")

#Load loops that filtered at the optimized cutoff (ipet.count >=6, FDR < 0.05).
#Analysis are performed only on the pooled samples, since the correlations between replicates are good.

loop.N <- cluster.N%>%filter(ipet.counts >=6 , p.adjust < 0.05, distance > 0)%>%
  select(header.loop)%>%mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))

loop.D <- cluster.D%>%filter(ipet.counts >=6 , p.adjust < 0.05, distance > 0)%>%
  select(header.loop)%>%mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))

loop.RE <- cluster.RE%>%filter(ipet.counts >=6 , p.adjust < 0.05, distance > 0)%>%
  select(header.loop)%>%mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))

loop.bzip23 <- cluster.bzip23%>%filter(ipet.counts >=6 , p.adjust < 0.05, distance > 0)%>%
  select(header.loop)%>%mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))

loop.N%>%filter(str_detect(loop.index,"Chr01_2254083_2280089")==T,type==0)%>%dim()
loop.D%>%filter(str_detect(loop.index,"Chr01_2254083_2280089")==T,type==0)%>%dim()
loop.RE%>%filter(str_detect(loop.index,"Chr01_2254083_2280089")==T,type==0)%>%dim()

####=========================================================================================================######
#Plot the cluster span distributions and pick appropriate length thresholds for view and CID analysis
####=========================================================================================================######

loop.length <- rbind(loop.N%>%filter(chrom1==chrom2)%>%select(distance)%>%mutate(condition="NC"),
                     loop.D%>%filter(chrom1==chrom2)%>%select(distance)%>%mutate(condition="DS"),
                     loop.RE%>%filter(chrom1==chrom2)%>%select(distance)%>%mutate(condition="RE"),
                     loop.bzip23%>%filter(chrom1==chrom2)%>%select(distance)%>%mutate(condition="bzip23")
)


loop.length$condition <- factor(loop.length$condition,levels = rev(c("NC","DS","RE","bzip23")),ordered = T)

#Anova
loop.length%>%mutate(log2.distance=log2(distance))%>%
  lm(log2.distance~condition,data = .)%>%
  aov()%>%
  TukeyHSD()

loop.length.hsd <- loop.length%>%mutate(log2.distance=log2(distance))%>%
  lm(log2.distance~condition,data = .)%>%
  glht(.,linfct = mcp(condition = "Tukey"))%>%cld(decreasing=T)

loop.length.hsd<-loop.length.hsd$mcletters$Letters%>%
  as.data.frame()%>%
  rename(hsd=".")%>%
  mutate(condition=rownames(.))

#Plot loop span distributions
pdf(file = "../04_loop_profiles/01_IntraChrom_loopSpan.distribute.pdf",
    height=2.28,width=3.72,bg=bg)

ggplot(loop.length,aes(x=log10(distance),
                       y=after_stat(count),
                       color=condition,group=condition))+
  geom_density(size=0.75)  +
  scale_color_manual(values = condition.col)+
  #facet_wrap(condition~.,nrow = 1,ncol = 3,scales = "fixed")+
  xlab("Log10(loop span)")+
  ylab("Loop count")

dev.off()


#Plot loop span distributions in boxplot
pdf(file = "../04_loop_profiles/02_IntraChrom_loopSpan.boxplot.pdf",
    height=2.28,width=3.72,bg=bg)

ggplot(loop.length,aes(y=log10(distance),x=condition))+
  geom_violin(aes(fill=condition),trim = F,scale="count",width=1)+
  geom_boxplot(fill="white",outlier.shape = NA,width=0.25)+
  stat_summary(fun.y = mean,size=0.1)+
  scale_fill_manual(values = condition.col)+
  geom_text(data=loop.length.hsd,
            aes(x=condition,y=1.05*log10(max(loop.length$distance)),label=hsd),
            inherit.aes = F)+
  ylab("Log10(loop span)") +
  coord_flip()

dev.off()


#Data summary
loop.length%>%group_by(condition)%>%summarize(median(log10(distance)))
loop.length%>%group_by(condition)%>%summarize(mean(log10(distance)))
loop.length%>%group_by(condition)%>%summarize(quantile(log10(distance)))


###Set loop length threshold for local loops:
#it is also the threshold for visible loops in tracks.
loopVizTH.D <- 200000
loopVizTH.N <- 200000
loopVizTH.RE <- 200000
loopVizTH.bzip23 <- 200000
####=========================================================================================================######
#Extract local/distant/interChr/interChr/all loops and output BEDPE files.
####=========================================================================================================######

for (x in sample.tag) {
  clust <- get(paste0("loop.",x))
  th <- get(paste0("loopVizTH.",x))

  #local loops with length < 200,000 bp
  write.table(file = paste0("../01_loop_BEDPE/",x,".local.loop.bedpe"),
              clust%>%dplyr::filter(distance <= th,type==1 )%>%
                dplyr::select(chrom1 ,start1 , end1 ,chrom2 ,start2  , end2,loop.index,ipet.counts),
              quote = F,row.names = F,col.names = F,sep = "\t")

  #distant loops with length >= 200,000 bp
  write.table(file = paste0("../01_loop_BEDPE/",x,".distant.loop.bedpe"),
              clust%>%dplyr::filter(distance > th,type==1 )%>%
                dplyr::select(chrom1 ,start1 , end1 ,chrom2 ,start2  , end2,loop.index,ipet.counts),
              quote = F,row.names = F,col.names = F,sep = "\t")
  #
  write.table(file = paste0("../01_loop_BEDPE/",x,".intraChr.loop.bedpe"),
              clust%>%dplyr::filter(type==1)%>%
                dplyr::select(chrom1 ,start1 , end1 ,chrom2 ,start2  , end2,loop.index,ipet.counts),
              quote = F,row.names = F,col.names = F,sep = "\t")

  write.table(file = paste0("../01_loop_BEDPE/",x,".interChr.loop.bedpe"),
              clust%>%dplyr::filter(type==0)%>%
                dplyr::select(chrom1 ,start1 , end1 ,chrom2 ,start2  , end2,loop.index,ipet.counts),
              quote = F,row.names = F,col.names = F,sep = "\t")

  write.table(file = paste0("../01_loop_BEDPE/",x,".all.loop.bedpe"),
              clust%>%dplyr::select(chrom1 ,start1 , end1 ,chrom2 ,start2  , end2,loop.index,ipet.counts),
              quote = F,row.names = F,col.names = F,sep = "\t")

  rm(clust,x)
}

####=========================================================================================================######
####annotate the associated genomic features (in respect to the nearest genes) of the clusters
####=========================================================================================================######

#Genes and features annotated in this section are used as materials for further analysis, not as direct results.
#The annotations can also be load from the ChIP-seq result directly.

library(GenomicFeatures)
library(ChIPseeker)
mh63.txdb <- makeTxDbFromGFF("~/genome_anno/MH63RS3/MH63RS3.gff3",format = "gff3")

annotateBEDPE <- function(x){
  anno <- mh63.encode%>%
    dplyr::select(MH63RS3_ID,MSU7_ID,name,Title,Description)%>%
    dplyr::rename(geneId=MH63RS3_ID)

  df <- x
  anchor1 <- df%>%
    dplyr::select(chrom1,start1 ,end1,loop.index)%>%
    makeGRangesFromDataFrame(seqnames.field = "chrom1",
                             start.field="start1",
                             end.field="end1",
                             keep.extra.columns=T)%>%
    annotatePeak(peak = .,tssRegion = c(-1000,+100),
                 TxDb = mh63.txdb,level = "gene",
                 flankDistance = 1000,
                 addFlankGeneInfo = T)%>%
    #as.data.frame()%>%dplyr::select(-seqnames,-start,-end,-strand)%>%
    as.data.frame()%>%
    dplyr::rename(anchor1.chr=seqnames,
                  anchor1.str=start,
                  anchor1.end=end)%>%
    dplyr::select(-strand)%>%
    merge(.,anno,by="geneId",all.x=T)%>%
    dplyr::rename(width1=width,
                  annotation1=annotation,
                  geneId1=geneId,
                  distanceToTSS1=distanceToTSS,
                  flank_geneIds1=flank_geneIds,
                  flank_gene_distances1=flank_gene_distances,
                  MSU7_ID1=MSU7_ID,
                  Gene_Symbol1=name,
                  Title1=Title,
                  Description1=Description,
                  geneChr1=geneChr,
                  geneStart1=geneStart,
                  geneEnd1=geneEnd,
                  geneStrand1=geneStrand,
                  geneLength1=geneLength)

  anchor2 <- df%>%
    dplyr::select(chrom2,start2 ,end2,loop.index,ipet.counts, type, distance)%>%
    makeGRangesFromDataFrame(seqnames.field = "chrom2",start.field="start2",end.field="end2",keep.extra.columns=T)%>%
    annotatePeak(peak = .,tssRegion = c(-1000,+100),
                 TxDb = mh63.txdb,level = "gene",
                 flankDistance = 1000,
                 addFlankGeneInfo = T)%>%
    #as.data.frame()%>%dplyr::select(-seqnames,-start,-end,-strand)%>%
    as.data.frame()%>%
    dplyr::rename(anchor2.chr=seqnames,anchor2.str=start,anchor2.end=end)%>%
    dplyr::select(-strand)%>%
    merge(.,anno,by="geneId",all.x=T)%>%
    dplyr::rename(width2=width,
                  annotation2=annotation,
                  geneId2=geneId,
                  distanceToTSS2=distanceToTSS,
                  flank_geneIds2=flank_geneIds,
                  flank_gene_distances2=flank_gene_distances,
                  MSU7_ID2=MSU7_ID,
                  Gene_Symbol2=name,
                  Title2=Title,
                  Description2=Description,
                  geneChr2=geneChr,
                  geneStart2=geneStart,
                  geneEnd2=geneEnd,
                  geneStrand2=geneStrand,
                  geneLength2=geneLength)

  full_join(anchor1,anchor2,by="loop.index")
}

for (x in sample.tag) {
  get(paste0("loop.",x))%>%
    annotateBEDPE(x=.)%>%assign(paste0("chipseeker.loop.",x),.,envir = .GlobalEnv)

  get(paste0("chipseeker.loop.",x))%>%
    write.xlsx(file = paste0("../02_Loop_feature_annotation/","all.loops.",x,".ChIPseeker.xlsx"),.,overwrite = T)

  rm(x)
}

detach(name = "package:GenomicFeatures")
detach(name = "package:AnnotationDbi") #detach AnnotationDbi to avoid conflict with dplyr

####=========================================================================================================######
#plot the proportions of each interacting feature (iFeature) combinations
####=========================================================================================================######

getPairedFeatures <- function(x){

  out <- x%>%mutate(

    iFeature.tmp = str_c(annotation1,annotation2),

    iFeature=ifelse(
      grepl(x=annotation1,pattern = "^Promoter")==T &
        grepl(x=annotation2,pattern = "^Promoter")==T,
      "Promoter-Promoter",
      ifelse(
        grepl(x=iFeature.tmp,pattern = "Promoter")==T &
          grepl(x=iFeature.tmp,pattern = "Distal")==T,
        "Promoter-Distal intergenic",
        ifelse(grepl(x=iFeature.tmp,pattern = "Promoter")==T &
                 (grepl(x=iFeature.tmp,pattern = "UTR")==T |
                    grepl(x=iFeature.tmp,pattern = "Exon")==T |
                    grepl(x=iFeature.tmp,pattern = "Intron")==T ),
               "Promoter-Gene body",
               "Other"
        )
      )
    )
  )%>%
    select(loop.index,type,distance,annotation1,annotation2,iFeature )
  return(out)
}

iFeature.loop.H3K9ac.D <- getPairedFeatures(chipseeker.loop.D)
iFeature.loop.H3K9ac.N <- getPairedFeatures(chipseeker.loop.N)
iFeature.loop.H3K9ac.RE <- getPairedFeatures(chipseeker.loop.RE)
iFeature.loop.bzip23 <- getPairedFeatures(chipseeker.loop.bzip23)


#Plot iFeature proportion
pdf(file = "../04_loop_profiles/03_iFeature_proportion.pdf",height=2,width=3.72,bg=bg)

rbind(iFeature.loop.H3K9ac.D%>%mutate(condition="DS"),
      iFeature.loop.H3K9ac.N%>%mutate(condition="NC"),
      iFeature.loop.H3K9ac.RE%>%mutate(condition="RE"),
      iFeature.loop.bzip23%>%mutate(condition="bzip23"))%>%
  mutate(condition=factor(condition,levels = c("NC","DS","RE","bzip23"),ordered = T))%>%
  ggplot(aes(x=condition,fill=iFeature))+
  geom_histogram(position = "fill",stat = "count")+
  scale_fill_manual(values = c("Promoter-Promoter"="#7e0100",
                               "Promoter-Gene body"="#fea318",
                               "Promoter-Distal intergenic"="#898e43",
                               "Other"="#145e82"))+
  ylab("Percentage loops")+
  xlab("")

dev.off()

####=========================================================================================================######
#### Unfold the flanking genes in the annotation object to get full list of PPI genes
####=========================================================================================================######
#Promoter-featured anchors often overlap with multiple genes or associated with multiple adjacent promoters.
#In case of the need for analysis on all anchor-covered/associated genes, we in this section unfold the flanking
#genes within the promoter range annotated by ChipSeeker.
#This conversion only inculdes promoter-promoter interacting loops, genes associated to other loops are discarded.

#Load WGCNA modules classified based on RiceXpro seedling phytohormone treatment microarray.
wgcna.mod <- read.xlsx("../../03_RNAseq_MH63/07_Gene_categories_by_expression/all_gene_categories.xlsx")%>%
  select(-Dynamic.to.DS)

#Load mfuzz co-expression modules by MH63-WT RNA-seq.
#Module 0 means no module is assigned to the gene.
mfuzz.clust <- read.xlsx("../../03_RNAseq_MH63/05_CoExpression_network/02_mfuzz/ExpressedGenes.mfuzzModuleMemship.xlsx")%>%
  rename(geneID=gid)%>%
  select(-Row.names)%>%
  rename(membership.mfuzz.1="1",
         membership.mfuzz.2="2",
         membership.mfuzz.3="3",
         membership.mfuzz.4="4",
         membership.mfuzz.5="5",
         membership.mfuzz.6="6")

gene.cat <- merge(wgcna.mod,mfuzz.clust,by="geneID",all.x=T)

gene.cat[is.na(gene.cat)] <- 0

#Function to extract all promoter-associated or covered genes by anchors and add annotations.
get.PPI.gene.all <- function(x){
  anno <- mh63.encode%>%dplyr::select(MH63RS3_ID,MSU7_ID,name,Keyword,Description,Title)

  anchor1 <- x%>%filter(annotation1=="Promoter",annotation2=="Promoter")%>%
    select(loop.index,anchor1.chr,anchor1.str,anchor1.end,flank_geneIds1,flank_gene_distances1)%>%
    tidyr::separate_rows(flank_geneIds1,flank_gene_distances1,sep=";",convert=T)%>%as.data.frame%>%
    filter(flank_gene_distances1 > -1000, flank_gene_distances1 <= 0)%>%unique()%>%
    merge(.,anno,by.x="flank_geneIds1",by.y="MH63RS3_ID",all.x=T)%>%
    merge(.,gene.cat,by.x="flank_geneIds1",by.y="geneID",all.x=T)%>%
    mutate(anchor1.id=str_c(anchor1.chr,anchor1.str,anchor1.end,sep = "_"))%>%
    dplyr::rename(MSU7_ID.gene1=MSU7_ID,Symbol.gene1=name,Keyword.gene1=Keyword,
                  Description.gene1=Description,Title.gene1=Title,wgcna.mod.gene1=CEM,
                  mfuzz.clust.gene1=mfuzz.cluster,expressed.gene1 = exped,
                  DEG.DS_vs_NC.gene1 =DEG.DS_vs_NC,DEG.RE_vs_DS.gene1=DEG.RE_vs_DS)

  anchor2 <- x%>%filter(annotation1=="Promoter",annotation2=="Promoter")%>%
    select(loop.index,anchor2.chr,anchor2.str,anchor2.end,flank_geneIds2,flank_gene_distances2,ipet.counts,type)%>%
    tidyr::separate_rows(flank_geneIds2,flank_gene_distances2,sep=";",convert=T)%>%as.data.frame%>%
    filter(flank_gene_distances2 > -1000, flank_gene_distances2 <= 0)%>%unique()%>%
    merge(.,anno,by.x="flank_geneIds2",by.y="MH63RS3_ID",all.x=T)%>%
    merge(.,gene.cat,by.x="flank_geneIds2",by.y="geneID",all.x=T)%>%
    mutate(anchor2.id=str_c(anchor2.chr,anchor2.str,anchor2.end,sep = "_"))%>%
    dplyr::rename(MSU7_ID.gene2=MSU7_ID,Symbol.gene2=name,Keyword.gene2=Keyword,
                  Description.gene2=Description,Title.gene2=Title,wgcna.mod.gene2=CEM,
                  mfuzz.clust.gene2=mfuzz.cluster,expressed.gene2 = exped,
                  DEG.DS_vs_NC.gene2 =DEG.DS_vs_NC,DEG.RE_vs_DS.gene2=DEG.RE_vs_DS)

  tmp <- merge(anchor1,anchor2,by="loop.index")%>%
    arrange(anchor1.chr,anchor2.chr,anchor1.str,anchor2.str,flank_geneIds1,flank_geneIds2)

  tmp$wgcna.mod.gene1[is.na(tmp$wgcna.mod.gene1)] <- "grey"

  tmp$wgcna.mod.gene2[is.na(tmp$wgcna.mod.gene2)] <- "grey"


  return(tmp)
}

#get the PPI loop annotations with all associated genes.
for (x in sample.tag) {
  get(paste0("chipseeker.loop.",x))%>%
    get.PPI.gene.all(x=.)%>%
    assign(paste0("PPIloop.allGenes.",x),.,envir = .GlobalEnv)

  get(paste0("PPIloop.allGenes.",x))%>%
    write.xlsx(file = paste0("../03_PPI_pairs/01_PPI_loops_with_all_genes/",x,".PPI.xlsx"),.,overwrite = T)

  rm(x)
}

#Output the information of PPI genes in all categories for network node annotation and facet plot.

rbind(
  PPIloop.allGenes.D%>%
    select(anchor1.id,flank_geneIds1,flank_gene_distances1)%>%
    rename(anchor.id=anchor1.id,flank_geneIds=flank_geneIds1,flank_gene_distances=flank_gene_distances1),

  PPIloop.allGenes.D%>%
    select(anchor2.id,flank_geneIds2,flank_gene_distances2)%>%
    rename(anchor.id=anchor2.id,flank_geneIds=flank_geneIds2,flank_gene_distances=flank_gene_distances2),

  PPIloop.allGenes.N%>%
    select(anchor1.id,flank_geneIds1,flank_gene_distances1)%>%
    rename(anchor.id=anchor1.id,flank_geneIds=flank_geneIds1,flank_gene_distances=flank_gene_distances1),

  PPIloop.allGenes.N%>%
    select(anchor2.id,flank_geneIds2,flank_gene_distances2)%>%
    rename(anchor.id=anchor2.id,flank_geneIds=flank_geneIds2,flank_gene_distances=flank_gene_distances2),

  PPIloop.allGenes.RE%>%
    select(anchor1.id,flank_geneIds1,flank_gene_distances1)%>%
    rename(anchor.id=anchor1.id,flank_geneIds=flank_geneIds1,flank_gene_distances=flank_gene_distances1),

  PPIloop.allGenes.RE%>%
    select(anchor2.id,flank_geneIds2,flank_gene_distances2)%>%
    rename(anchor.id=anchor2.id,flank_geneIds=flank_geneIds2,flank_gene_distances=flank_gene_distances2),

  PPIloop.allGenes.bzip23%>%
    select(anchor1.id,flank_geneIds1,flank_gene_distances1)%>%
    rename(anchor.id=anchor1.id,flank_geneIds=flank_geneIds1,flank_gene_distances=flank_gene_distances1),

  PPIloop.allGenes.bzip23%>%
    select(anchor2.id,flank_geneIds2,flank_gene_distances2)%>%
    rename(anchor.id=anchor2.id,flank_geneIds=flank_geneIds2,flank_gene_distances=flank_gene_distances2)

)%>%
  unique()%>%
  merge(.,mh63.encode%>%dplyr::select(MH63RS3_ID,MSU7_ID,name,Keyword,Description,Title),
        by.x="flank_geneIds",by.y="MH63RS3_ID",all.x=T)%>%
  merge(.,gene.cat,by.x="flank_geneIds",by.y="geneID",all.x=T)%>%
  merge(.,tpm.all,by.x="flank_geneIds",by.y="gene_id",all.x=T)%>%
  write.xlsx("../03_PPI_pairs/02_All_PPI_gene_anno/PPI_Gene.anno.allConditions.xlsx")

####=========================================================================================================######
####Compare the expression levels of inter-chrom and intra-chrom PPI genes.
####=========================================================================================================######

#Get the inter-chrom and intra-chrom gene lists
get.intraInter.PPIgenes <- function(x){

  intraChr.PPI.genes <- get(paste0("PPIloop.allGenes.",x))%>%
    rowwise%>%
    filter(str_split(anchor1.id,pattern="_",simplify = T)[1]==str_split(anchor2.id,pattern="_",simplify = T)[1])%>%
    select(flank_geneIds1,flank_geneIds2)

  intraChr.PPI.genes <- c(intraChr.PPI.genes$flank_geneIds1,intraChr.PPI.genes$flank_geneIds2)%>%unique()


  interChr.PPI.genes <- get(paste0("PPIloop.allGenes.",x))%>%
    rowwise%>%
    filter(str_split(anchor1.id,pattern="_",simplify = T)[1] != str_split(anchor2.id,pattern="_",simplify = T)[1])%>%
    select(flank_geneIds1,flank_geneIds2)

  interChr.PPI.genes<- c(interChr.PPI.genes$flank_geneIds1,interChr.PPI.genes$flank_geneIds2)%>%unique()

  assign(paste0("intraChr.PPI.genes.",x),intraChr.PPI.genes,envir = .GlobalEnv)
  assign(paste0("interChr.PPI.genes.",x),interChr.PPI.genes,envir = .GlobalEnv)
}

lapply(c("D","N","RE"), get.intraInter.PPIgenes)


#We can expect that some genes are associated to intra and inter-chrom loops simultaneously.
#So, we investigate their relationship by upset plot.

library(UpSetR)
library(tidyverse)

for (x in c("D","N","RE")) {

  intra <- get(paste0("intraChr.PPI.genes.",x))
  inter <- get(paste0("interChr.PPI.genes.",x))

  rbind(data.frame(ID=intra,condition="intraChr"),
        data.frame(ID=inter,condition="interChr")
  )%>%
    as_tibble()%>%
    mutate(occur=1)%>%
    pivot_wider(
      id_cols = ID,
      names_from = condition,
      values_from = occur,
      values_fill = list(occur = 0)
    ) %>%
    data.frame()%>%
    left_join(tpm.all,by=c("ID"="gene_id"))%>%
    assign(paste0("upset.interIntraPPIgenes.",x),.,envir = .GlobalEnv)
}

rm(inter,intra)


#plot the upset plot together with the TPM value of the corresponding set.
#===============================================================================================
#N
pdf("../04_loop_profiles/04_interIntra_PPI.NC.upset.pdf",height =6,width = 3,bg=bg)
upset(upset.interIntraPPIgenes.N%>%
        mutate(log2TPMp1=log2(TPM.NC+1)),
      sets = c("intraChr","interChr"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      boxplot.summary = c("log2TPMp1"))
#sets.bar.color =c(col.NC,col.DS,col.RE))
dev.off()

upset.interIntraPPIgenes.N.hsd <- upset.interIntraPPIgenes.N%>%
  mutate(type=as.factor(case_when(
    intraChr==1 & interChr==0 ~ "set1",
    intraChr==1 & interChr==1 ~ "set2",
    intraChr==0 & interChr==1 ~ "set3"
  )))%>%
  lm(log2(TPM.NC+1)~type,data = .)%>%
  glht(.,linfct = mcp(type = "Tukey"))%>%cld(decreasing=T)


#D
pdf("../04_loop_profiles/05_interIntra_PPI.DS.upset.pdf",height = 6,width = 3,bg=bg)
upset(upset.interIntraPPIgenes.D%>%
        mutate(log2TPMp1=log2(TPM.DS+1)),
      sets = c("intraChr","interChr"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      boxplot.summary = c("log2TPMp1"))
#sets.bar.color =c(col.NC,col.DS,col.RE))
dev.off()

upset.interIntraPPIgenes.D.hsd <- upset.interIntraPPIgenes.D%>%
  mutate(type=as.factor(case_when(
    intraChr==1 & interChr==0 ~ "set1",
    intraChr==1 & interChr==1 ~ "set2",
    intraChr==0 & interChr==1 ~ "set3"
  )))%>%
  lm(log2(TPM.DS+1)~type,data = .)%>%
  glht(.,linfct = mcp(type = "Tukey"))%>%cld(decreasing=T)

#RE
pdf("../04_loop_profiles/06_interIntra_PPI.RE.upset.pdf",height = 6,width = 3,bg=bg)
upset(upset.interIntraPPIgenes.RE%>%
        mutate(log2TPMp1=log2(TPM.RE+1)),
      sets = c("intraChr","interChr"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      boxplot.summary = c("log2TPMp1"))
#sets.bar.color =c(col.NC,col.DS,col.RE))
dev.off()

upset.interIntraPPIgenes.RE.hsd <- upset.interIntraPPIgenes.RE%>%
  mutate(type=as.factor(case_when(
    intraChr==1 & interChr==0 ~ "set1",
    intraChr==1 & interChr==1 ~ "set2",
    intraChr==0 & interChr==1 ~ "set3"
  )))%>%
  lm(log2(TPM.RE+1)~type,data = .)%>%
  glht(.,linfct = mcp(type = "Tukey"))%>%cld(decreasing=T)


##In general inter-chr PPI add a significant plus to the transcription level of the anchored genes.

####=========================================================================================================######
####Compare the expression levels of PPI genes, non-PPI genes, basal genes, and unmarked genes under each condition.
####=========================================================================================================######

#####Compare the PPI genes (use all genes) to basal H3K9ac-marked genes and unmarked genes under each condition
#####(do not include bzip23 data here).

####Get basal genes and promoter-distal genes (promoter side only).
#load the list of all promoter-associated or covered genes by IDR peaks under the three conditions.
IDRpeak.proGene.D <- read.xlsx("../../01_ChIPseq_MH63/08_Promoter_marked_genes/IDRpeaks.DS.AllPromoterMarkedGenes.xlsx")%>%
  select(geneId)%>%unique() #17769
IDRpeak.proGene.N <- read.xlsx("../../01_ChIPseq_MH63/08_Promoter_marked_genes/IDRpeaks.NC.AllPromoterMarkedGenes.xlsx")%>%
  select(geneId)%>%unique() #17237
IDRpeak.proGene.RE <- read.xlsx("../../01_ChIPseq_MH63/08_Promoter_marked_genes/IDRpeaks.RE.AllPromoterMarkedGenes.xlsx")%>%
  select(geneId)%>%unique() #17604

#We defined that if a H3K9ac-marked gene is not associated to a loop (either a PPI loop or other types), it is a basal gene.

lapply(c("D","N","RE"),
       function(x){

         #PPI <- get(paste0("PPIloop.allGenes.",x))%>%select(flank_geneIds1, flank_geneIds2)
         #PPI.genes <- unique(c(PPI$flank_geneIds1,PPI$flank_geneIds2))
         PPI.genes <- get(paste0("chipseeker.loop.",x))%>%
           filter((annotation1== "Promoter" & annotation2=="Promoter"))
         PPI.genes <- unique(c(PPI.genes$geneId1,PPI.genes$geneId1))

         ProDistal.genes <- get(paste0("chipseeker.loop.",x))%>%
           filter((annotation1== "Promoter" & annotation2=="Distal Intergenic")|
                    (annotation2== "Promoter" & annotation1=="Distal Intergenic"))%>%
           #mutate(geneId=case_when(annotation1== "Promoter" ~ geneId1,
           #                        annotation2== "Promoter" ~ geneId2,))%>%

           mutate(flank_geneIds=case_when(annotation1== "Promoter" ~ flank_geneIds1,
                                          annotation2== "Promoter" ~ flank_geneIds2),
                  flank_gene_distances=case_when(annotation1== "Promoter" ~ flank_gene_distances1,
                                                 annotation2== "Promoter" ~ flank_gene_distances2))%>%
           tidyr::separate_rows(flank_geneIds,flank_gene_distances,sep=";",convert=T)%>%as.data.frame%>%
           filter(flank_gene_distances > -1000, flank_gene_distances <= 0)%>%
           select(flank_geneIds)%>%
           rename(geneId=flank_geneIds)%>%
           unique()

         ProDistal.genes <- ProDistal.genes$geneId

         basalGenes <- get(paste0("IDRpeak.proGene.",x))%>%
           filter(!(geneId %in% c(PPI.genes,ProDistal.genes)))%>%
           unique()

         basalGenes <-  basalGenes$geneId

         unmarked <- tpm.all%>%
           filter(!(gene_id %in% c(PPI.genes,ProDistal.genes,basalGenes)))%>%
           unique()

         unmarked <- unmarked$gene_id

         assign(paste0("PPIgenes.",x),PPI.genes,envir = .GlobalEnv)
         assign(paste0("ProDistalGenes.",x),ProDistal.genes,envir = .GlobalEnv)
         assign(paste0("basalGenes.",x),basalGenes,envir = .GlobalEnv)
         assign(paste0("unmarkedGenes.",x),unmarked,envir = .GlobalEnv)

       }
)

loop.class.TPM.data <- rbind(
  #DS
  tpm.all%>%mutate(iclass=case_when(
    gene_id%in%PPIgenes.D ~ "PPI",
    gene_id%in%ProDistalGenes.D ~ "ProDistal",
    gene_id%in%basalGenes.D ~ "basal",
    gene_id%in%unmarkedGenes.D ~ "unmarked"
  ))%>%
    select(gene_id,iclass,TPM.DS)%>%
    mutate(condition="DS")%>%
    dplyr::rename(TPM=TPM.DS),
  #NC
  tpm.all%>%mutate(iclass=case_when(
    gene_id%in%PPIgenes.N ~ "PPI",
    gene_id%in%ProDistalGenes.N ~ "ProDistal",
    gene_id%in%basalGenes.N ~ "basal",
    gene_id%in%unmarkedGenes.N ~ "unmarked"
  ))%>%
    select(gene_id,iclass,TPM.NC)%>%
    mutate(condition="NC")%>%
    dplyr::rename(TPM=TPM.NC),
  #RE
  tpm.all%>%mutate(iclass=case_when(
    gene_id%in%PPIgenes.RE ~ "PPI",
    gene_id%in%ProDistalGenes.RE ~ "ProDistal",
    gene_id%in%basalGenes.RE ~ "basal",
    gene_id%in%unmarkedGenes.RE ~ "unmarked"
  ))%>%
    select(gene_id,iclass,TPM.RE)%>%
    mutate(condition="RE")%>%
    dplyr::rename(TPM=TPM.RE)
)%>%mutate(condition=factor(condition,levels = c("NC","DS","RE"),ordered = T))%>%
  mutate(iclass=factor(iclass,levels = c("PPI","ProDistal","basal","unmarked"),ordered = T))

pdf("../04_loop_profiles/07_PPI_basal_expression_compare.violin.pdf",width = 10.5,height = 4,bg=bg)
loop.class.TPM.data%>%
  ggplot(aes(x=iclass,y=log2(TPM+1)))+
  geom_violin(aes(fill=condition),trim = F,scale = "width")+
  geom_boxplot(width=0.3,outlier.shape = NA)+
  scale_fill_manual(values = condition.col)+
  facet_wrap(.~condition,ncol = 3,scales = "fixed")+
  ylab("log2(TPM+1)")+
  xlab("Interaction categories")
dev.off()

#HSD
loop.class.TPM.hsd <- loop.class.TPM.data%>%
  mutate(type=as.factor(str_c(iclass,condition,sep = "_")))%>%
  lm(log2(TPM+1)~type,data = .)%>%
  glht(.,linfct = mcp(type = "Tukey"))%>%
  cld(decreasing=T)
#Export the IDs of PPI, basal, and unmarked genes under each conditions

write.xlsx(loop.class.TPM.data,"../05_PPI_PDI_and_basal_genes/genes_by_loop_class.xlsx")


####=========================================================================================================######
#Plot the proportions of the expressed and non-expressed PPI genes.
####=========================================================================================================######

#In this analysis, all PPI genes are used.
pdf("../04_loop_profiles/08_PPI_gene_expressed_percentage.pdf",height = 3,width = 3,bg=bg)
rbind(
  data.frame(PPI.wgcna.mod=c(PPIloop.allGenes.D$expressed.gene1,PPIloop.allGenes.D$expressed.gene2),condition="DS"),
  data.frame(PPI.wgcna.mod=c(PPIloop.allGenes.N$expressed.gene1,PPIloop.allGenes.N$expressed.gene2),condition="NC"),
  data.frame(PPI.wgcna.mod=c(PPIloop.allGenes.RE$expressed.gene1,PPIloop.allGenes.RE$expressed.gene2),condition="RE")
)%>%ggplot(aes(x=condition,fill=PPI.wgcna.mod))+
  geom_histogram(position = "fill",stat = "count")+
  ylab("Percentage PPI genes")
dev.off()

pdf("../04_loop_profiles/09_PPI_loop_expression_class.pdf",height = 3,width = 3,bg=bg)
rbind(
  PPIloop.allGenes.D%>%
    select(loop.index,expressed.gene1,expressed.gene2)%>%
    mutate(pair.ident=ifelse(expressed.gene1=="No" & expressed.gene2=="No", "dual-none",
                             ifelse(expressed.gene1=="No" | expressed.gene2=="No", "single-none","dual-exp")))%>%
    mutate(condition="DS"),

  PPIloop.allGenes.N%>%
    select(loop.index,expressed.gene1,expressed.gene2)%>%
    mutate(pair.ident=ifelse(expressed.gene1=="No" & expressed.gene2=="No", "dual-none",
                             ifelse(expressed.gene1=="No" | expressed.gene2=="No", "single-none","dual-exp")))%>%
    mutate(condition="NC"),

  PPIloop.allGenes.RE%>%
    select(loop.index,expressed.gene1,expressed.gene2)%>%
    mutate(pair.ident=ifelse(expressed.gene1=="No" & expressed.gene2=="No", "dual-none",
                             ifelse(expressed.gene1=="No" | expressed.gene2=="No", "single-none","dual-exp")))%>%
    mutate(condition="RE")

)%>%

  ggplot(aes(fill=pair.ident,x=condition))+geom_histogram(stat = "count",position = "fill")+
  ylab("Percentage PPI loops")

dev.off()

####=========================================================================================================######
####Plot inter-chrom interactions in circos layout.
####=========================================================================================================######

library(circlize)

pdf("../04_loop_profiles/10_DS_interChromLoops.circos.pdf",width =5,height = 5,bg =bg)
circos.initializeWithIdeogram(cytoband = "~/genome_anno/MH63RS3/HM63RS3.cytoband.txt" )
circos.genomicLink(loop.D%>%filter(chrom1 != chrom2)%>%select(chrom1,start1,end1),
                   loop.D%>%filter(chrom1 != chrom2)%>%select(chrom2,start2,end2),
                   col = alpha(col.DS,alpha  = 0.2),lwd=0.01)
circos.clear()
dev.off()

pdf("../04_loop_profiles/11_NC_interChromLoops.circos.pdf",width =5,height = 5,bg =bg)
circos.initializeWithIdeogram(cytoband = "~/genome_anno/MH63RS3/HM63RS3.cytoband.txt" )
circos.genomicLink(loop.N%>%filter(chrom1 != chrom2)%>%select(chrom1,start1,end1),
                   loop.N%>%filter(chrom1 != chrom2)%>%select(chrom2,start2,end2),
                   col = alpha(col.NC,alpha = 0.2),lwd=0.01)
#highlight the FISH targets
#circos.genomicLink(loop.N%>%filter(loop.index== "Chr01_2254083_2280089_Chr05_25484182_25501247")%>%select(chrom1,start1,end1),
#                   loop.N%>%filter(loop.index== "Chr01_2254083_2280089_Chr05_25484182_25501247")%>%select(chrom2,start2,end2),
#                   col = "magenta",lwd=1)

#circos.genomicLink(loop.H3K9ac.NC%>%filter(loop.index== "Chr02_24120421_24134330_Chr04_33033504_33046659")%>%select(chrom1,start1,end1),
#                   loop.H3K9ac.NC%>%filter(loop.index== "Chr02_24120421_24134330_Chr04_33033504_33046659")%>%select(chrom2,start2,end2),
#                   col = "magenta",lwd=0.1)

#circos.genomicLink(DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V1,V2,V3),
#                   DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V4,V5,V6),
#                   col = "magenta",lwd=1)

circos.clear()
dev.off()


pdf("../04_loop_profiles/12_RE_interChromLoops.circos.pdf",width =5,height = 5,bg =bg)
circos.initializeWithIdeogram(cytoband = "~/genome_anno/MH63RS3/HM63RS3.cytoband.txt" )
circos.genomicLink(loop.RE%>%filter(chrom1 != chrom2)%>%select(chrom1,start1,end1),
                   loop.RE%>%filter(chrom1 != chrom2)%>%select(chrom2,start2,end2),
                   col = alpha(col.RE,alpha = 0.2),lwd=0.01)

#circos.genomicLink(DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V1,V2,V3),
#                   DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V4,V5,V6),
#                   col = "magenta",lwd=1)

circos.clear()
dev.off()

pdf("../04_loop_profiles/13_bzip23_interChromLoops.circos.pdf",width =5,height = 5,bg =bg)
circos.initializeWithIdeogram(cytoband = "~/genome_anno/MH63RS3/HM63RS3.cytoband.txt" )
circos.genomicLink(loop.bzip23%>%filter(chrom1 != chrom2)%>%select(chrom1,start1,end1),
                   loop.bzip23%>%filter(chrom1 != chrom2)%>%select(chrom2,start2,end2),
                   col = alpha(col.bzip23,alpha = 0.2),lwd=0.01)

#circos.genomicLink(DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V1,V2,V3),
#                   DS.domPPIloop.bedpe%>%filter(V7 == "loop_DS_4803")%>%select(V4,V5,V6),
#                   col = "magenta",lwd=1)

circos.clear()
dev.off()

####=========================================================================================================######
####Plot inter-chrom interactions in synteny style.
####=========================================================================================================######
library(RIdeogram)

#Generate ideogram file
mh63.synteny.ideogram <- rbind(mh63.chrlength%>%
                                 mutate(Start=1,
                                        fill=252525,
                                        species="mh63.1",
                                        size=10,
                                        color=252525),
                               mh63.chrlength%>%
                                 mutate(Start=1,
                                        fill=252525,
                                        species="mh63.2",
                                        size=10,
                                        color=252525)
)%>%select(Chr,Start,End,fill,species,size,color)

lapply(c("N","D","RE","bzip23"), function(x){

  get(paste0("loop.",x))%>%
    filter(chrom1 != chrom2)%>%
    rowwise()%>%
    mutate(index1=str_c(chrom1,start1,end1,sep = "_"),
           index2=str_c(chrom2,start2,end2,sep = "_"))%>%
    mutate(chrom1=str_replace(chrom1,pattern = "Chr",replacement = ""),
           chrom2=str_replace(chrom2,pattern = "Chr",replacement = ""))%>%
    mutate(chrom1=as.numeric(chrom1),chrom2=as.numeric(chrom2))%>%   #Chr names in the file should be numeric.
    mutate(index.a1 = case_when(chrom1<chrom2 ~ index1, chrom1>chrom2 ~ index2),
           index.a2 = case_when(chrom1<chrom2 ~ index2, chrom1>chrom2 ~ index1))%>%
    mutate(Species_1=str_split(index.a1,pattern = "_",simplify = T)[1],
           Start_1=str_split(index.a1,pattern = "_",simplify = T)[2]%>%as.numeric(),
           End_1=str_split(index.a1,pattern = "_",simplify = T)[3]%>%as.numeric(),
           Species_2=str_split(index.a2,pattern = "_",simplify = T)[1],
           Start_2=str_split(index.a2,pattern = "_",simplify = T)[2]%>%as.numeric(),
           End_2=str_split(index.a2,pattern = "_",simplify = T)[3]%>%as.numeric())%>%
    select(Species_1,Start_1,End_1,Species_2,Start_2,End_2)%>%
    mutate(Species_1=str_replace(Species_1,pattern = "Chr",replacement = ""),
           Species_2=str_replace(Species_2,pattern = "Chr",replacement = ""))%>%
    mutate(Species_1=as.numeric(Species_1),Species_2=as.numeric(Species_2))%>%
    mutate(fill="cccccc")%>%
    as.data.frame()%>%
    arrange(Species_1,Species_2,Start_1,Start_2)%>%
    assign(paste0("interChr.ideogram.",x),.,envir = .GlobalEnv)

})

#plot the interChr interactions between Chr01 and Chr05 and highlight the FISH sites (Chr01:2254083-2280089 // Chr05:25484182-25501247)


lapply(c("N","D","RE","bzip23"), function(x){

  ka <- rbind(mh63.synteny.ideogram%>%filter((Chr=="Chr01" & species=="mh63.1")),
              data.frame(Chr="pseudo1",     #make two pseudo intervals to occupy the flanking spaces.
                         Start=1,
                         End=6859802,
                         fill=252525,
                         species="mh63.2",
                         size=10,
                         color=252525),
              mh63.synteny.ideogram%>%filter((Chr=="Chr05" & species=="mh63.2")),
              data.frame(Chr="pseudo2",
                         Start=1,
                         End=6859802,
                         fill=252525,
                         species="mh63.2",
                         size=10,
                         color=252525)
  )


  sy <- get(paste0("interChr.ideogram.",x))%>%
    filter(Species_1 ==1,Species_2 ==5)%>%
    mutate(Species_2 =2)%>% #The real Chr05 is in the center (site 2).
    mutate(fill=case_when(Start_1==2254083 & Start_2 ==25484182 ~ "fcb06b",
                          TRUE ~ "cccccc"))
  ou <- paste0("../04_loop_profiles/14_",x,"_InterChrmLoops.synteny.svg")

  ideogram(karyotype = ka,
           synteny = sy,
           output= ou)

})

#The output should be carefully adjuested in AI.

####=========================================================================================================######
####The representative gene of anchors and PPI network analysis.
####=========================================================================================================######
#For clear network labeling, we use the anchor-associated gene with the highest (average) expression to represent
#the anchor. This can avoid error in node labeling if the anchor covers multiple genes in different categories.

#Generate the list for the unified anchor. The average TPM of all genes associated to the anchor is also preserved
#in case the TPM of the anchor should be labeled.
#For category information labeling of the node (anchor), the information of the representative gene is used.

#load DEG fold changes for annotation.

DE.DS_vs_NC <- read.xlsx("../../03_RNAseq_MH63/02_DEG/02_DEG_spreadsheet/allGenes_MH_D_vs_MH_N.xlsx")%>%
  select("Row.names","log2FoldChange")%>%rename(flank_geneIds=Row.names,log2FC_DS_vs_NC=log2FoldChange)

DE.RE_vs_DS <- read.xlsx("../../03_RNAseq_MH63/02_DEG/02_DEG_spreadsheet/allGenes_MH_RE_vs_MH_D.xlsx")%>%
  select("Row.names","log2FoldChange")%>%rename(flank_geneIds=Row.names,log2FC_RE_vs_DS=log2FoldChange)

DE.RE_vs_NC <- read.xlsx("../../03_RNAseq_MH63/02_DEG/02_DEG_spreadsheet/allGenes_MH_RE_vs_MH_N.xlsx")%>%
  select("Row.names","log2FoldChange")%>%rename(flank_geneIds=Row.names,log2FC_RE_vs_NC=log2FoldChange)

DE.DS_vs_bzip23 <- read.xlsx("../../26_DS_RNAseq_MH63_vs_bzip23/02_DEG/02_DEG_spreadsheet/allGenes_MH_DS_vs_bzip23_DS.xlsx")%>%
  select("Row.names","log2FoldChange")%>%rename(flank_geneIds=Row.names,log2FC_DS_vs_bzip23=log2FoldChange)


DE.total <- left_join(DE.DS_vs_NC,DE.RE_vs_DS)%>%left_join(.,DE.RE_vs_NC)%>%left_join(.,DE.DS_vs_bzip23)

#genes with non-quantified log2FCs are set to 0 to avoid error in averaging.

DE.total[is.na(DE.total)] <- 0

###annotate all anchors in all samples with the unique representative gene and the corresponding info.
#As all loops were called on a unified set of anchors, data from different samples can be r-bind directly
#to form a anchor info object.

PPI.anchor.info <- rbind(PPIloop.allGenes.D,
                         PPIloop.allGenes.N,
                         PPIloop.allGenes.RE,
                         PPIloop.allGenes.bzip23)%>%
  select(anchor1.id,anchor2.id,flank_geneIds1,flank_geneIds2)%>%
  unique()

PPI.anchor.info <- rbind(PPI.anchor.info%>%select(anchor1.id,flank_geneIds1)%>%
                           rename(anchor.id=anchor1.id,flank_geneIds=flank_geneIds1),
                         PPI.anchor.info%>%select(anchor2.id,flank_geneIds2)%>%
                           rename(anchor.id=anchor2.id,flank_geneIds=flank_geneIds2))%>%
  unique()

PPI.anchor.info %<>% left_join(.,tpm.all,by=c("flank_geneIds"="gene_id"))
PPI.anchor.info %<>% left_join(.,DE.total,by=c("flank_geneIds"))

PPI.anchor.info[is.na(PPI.anchor.info)] <- 0 #Non-quantified FC values are regarded as unchanged.

#filter by the max TPM to get the most active gene associated to the anchor,
#as H3K9ac is a typical activate mark, the most active gene can well represent the anchor.
#bzip23.TPM is not used in the filtering.

PPI.anchor.info.uniq <-  PPI.anchor.info%>%rowwise()%>%
  mutate(max.TPM.of.three=max(TPM.DS,TPM.NC,TPM.RE))%>%
  group_by(anchor.id)%>% #Next operations are within each anchor.
  mutate(DS.TPM.anchor=mean(TPM.DS), #Use the average TPM of all associated genes as the anchor TPM.
         NC.TPM.anchor=mean(TPM.NC),
         RE.TPM.anchor=mean(TPM.RE),
         bzip23.TPM.anchor=mean(TPM.bzip23),
         log2FC_DS_vs_NC.anchor=mean(log2FC_DS_vs_NC),
         log2FC_RE_vs_DS.anchor=mean(log2FC_RE_vs_DS),
         log2FC_RE_vs_NC.anchor=mean(log2FC_RE_vs_NC),
         log2FC_DS_vs_bzip23.anchor=mean(log2FC_DS_vs_bzip23)
  )%>%
  filter(max.TPM.of.three==max(max.TPM.of.three))%>% #get the most active gene
  as.data.frame()%>%unique()

#deal with the duplicated rows with none expressed genes (the max.TPM are equal)
#select a random genes to represent the duplicated rows, since the none expressed genes do not make much difference.
PPI.anchor.info.uniq.u <- filter(PPI.anchor.info.uniq,!(anchor.id %in%  PPI.anchor.info.uniq$anchor.id[duplicated(PPI.anchor.info.uniq$anchor.id)]))
PPI.anchor.info.uniq.d <- filter(PPI.anchor.info.uniq,anchor.id %in%  PPI.anchor.info.uniq$anchor.id[duplicated(PPI.anchor.info.uniq$anchor.id)])%>%
  group_by(anchor.id)%>%
  filter(flank_geneIds==min(flank_geneIds))
#check size
PPI.anchor.info$anchor.id%>%unique()%>%length() == dim(PPI.anchor.info.uniq.u)[1] + dim(PPI.anchor.info.uniq.d)[1]

#The final anchor info with unique representative gene.
PPI.anchor.info.uniq <- rbind(PPI.anchor.info.uniq.u,PPI.anchor.info.uniq.d)%>%
  arrange(flank_geneIds)%>%
  select(-max.TPM.of.three,-TPM.DS, -TPM.NC,-TPM.RE, -TPM.bzip23,
         -log2FC_DS_vs_NC, -log2FC_RE_vs_DS, -log2FC_RE_vs_NC, -log2FC_DS_vs_bzip23)%>%
  rename(represent.gene = flank_geneIds)

rm(PPI.anchor.info.uniq.u,PPI.anchor.info.uniq.d)

#add co-expression clusters
PPI.anchor.info.uniq %<>% left_join(.,mfuzz.clust,by=c("represent.gene"="geneID"))%>%
  left_join(.,wgcna.mod,by=c("represent.gene"="geneID"))%>%
  rename(mfuzz.cluster.anchor=mfuzz.cluster,
         expressed.anchor=exped,
         DEG.DS_vs_NC.anchor=DEG.DS_vs_NC,
         DEG.RE_vs_DS.anchor=DEG.RE_vs_DS,
         WGCNA.anchor=CEM)

#Combined the information lists.
PPI.anchor.info %<>% left_join(.,PPI.anchor.info.uniq)

#output the anchor list
write.xlsx(PPI.anchor.info,"../03_PPI_pairs/03_PPI_anchor_annoation/PPI_anchor_wi_all_genes.annoation.xlsx")
write.xlsx(PPI.anchor.info.uniq,"../03_PPI_pairs/03_PPI_anchor_annoation/PPI_anchor_unique_genes.annoation.xlsx")

#annotate PPIs with representative genes and info.
lapply(c("N","D","RE","bzip23"),function(x){

  name.org <- colnames(PPI.anchor.info.uniq)[-1]
  name.1 <- str_c(name.org,".1")
  name.2 <- str_c(name.org,".2")

  anno.org <- colnames(mh63.encode)[-1]
  anno.1 <- str_c(anno.org,".1")
  anno.2 <- str_c(anno.org,".2")

  get(paste0("PPIloop.allGenes.",x))%>%
    select(loop.index,ipet.counts,type,
           anchor1.chr,anchor1.str,anchor1.end,anchor1.id,
           anchor2.chr,anchor2.str,anchor2.end,anchor2.id,
    )%>%unique()%>%
    #add TPM under each conditions to flanking genes.
    left_join(PPI.anchor.info.uniq,by=c("anchor1.id"="anchor.id"))%>%
    rename_with( ~ name.1,all_of(name.org))%>%
    left_join(mh63.encode,by=c("represent.gene.1"="MH63RS3_ID"))%>%
    rename_with( ~ anno.1,all_of(anno.org))%>%
    left_join(PPI.anchor.info.uniq,by=c("anchor2.id"="anchor.id"))%>%
    rename_with( ~ name.2,all_of(name.org))%>%
    left_join(mh63.encode,by=c("represent.gene.2"="MH63RS3_ID"))%>%
    rename_with( ~ anno.2,all_of(anno.org))%>%
    assign(paste0("PPIloops.uniqGene.",x),.,envir = .GlobalEnv)

  write.xlsx(get(paste0("PPIloops.uniqGene.",x)),
             paste0("../03_PPI_pairs/04_PPI_loops_with_representative_genes/",x,
                    ".PPI.uniqGene.xlsx"))

})


####=========================================================================================================######
####Plot GW view for intra-Chrom PPI loops.
####=========================================================================================================######
#The plotting script is modified from the LGL script in ChIA-PET tools V3 pipeline.
#load plot script
library(grid)
source("IntraChromPlotter.R")

#Cytoband file labeling chromosome features (blank regions are labeled as "gneg",
#centromeres are labled as "acen", other features are labeled as "gpos100").
cyto <- read.table("~/genome_anno/MH63RS3/HM63RS3.cytoband.txt", fill = T, header = F)

#Genomic features to be ploted as histograms. A 4-column file with the value in the 4th column
#Use the fold-enrichment by macs2 (pooled samples) as values.
#Load peaks from the pooled narrowPeaks and extract the regions that appeared in the IDR result
#Some narrowPeaks may have multiple values (due to multiple peak summits in the region), leave it alone anyway.
library(GenomicRanges)
#DS
peaks.DS.IDR <- readPeakFile("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_D_H3K9ac.IDR.uniq.narrowPeak")

peaks.DS <- read.table("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/01_narrowPeak/MH_D_H3K9ac_pooled_peaks.narrowPeak",sep = "\t",header = F)%>%
  dplyr::select(V1,V2,V3,V7)%>%unique()%>%mutate(V7=log2(V7))%>%makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)

peaks.DS <- peaks.DS[findOverlaps(peaks.DS,peaks.DS.IDR)@from,]%>%as.data.frame()%>%select(seqnames,start,end,V7)%>%filter(V7 > 1)

#NC
peaks.NC.IDR <- readPeakFile("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_N_H3K9ac.IDR.uniq.narrowPeak")

peaks.NC <- read.table("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/01_narrowPeak/MH_N_H3K9ac_pooled_peaks.narrowPeak",sep = "\t",header = F)%>%
  dplyr::select(V1,V2,V3,V7)%>%unique()%>%mutate(V7=log2(V7))%>%makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)

peaks.NC <- peaks.NC[findOverlaps(peaks.NC,peaks.NC.IDR)@from,]%>%as.data.frame()%>%select(seqnames,start,end,V7)%>%filter(V7 > 1)

#RE
peaks.RE.IDR <- readPeakFile("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/03_IDR/02_IDR_nonredundent_regions/MH_RE_H3K9ac.IDR.uniq.narrowPeak")

peaks.RE <- read.table("../../../01_ChIPseq/01_MH_ChIPseq/01_Peak_set/01_narrowPeak/MH_RE_H3K9ac_pooled_peaks.narrowPeak",sep = "\t",header = F)%>%
  dplyr::select(V1,V2,V3,V7)%>%unique()%>%mutate(V7=log2(V7))%>%makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)

peaks.RE <- peaks.RE[findOverlaps(peaks.RE,peaks.RE.IDR)@from,]%>%as.data.frame()%>%select(seqnames,start,end,V7)%>%filter(V7 > 1)

#Normalize the peak value to 0-1. Now the value means the histogram height relative to a y-axis ranging from the min log2 enrichment to the max.
peaks.DS[[4]] <- (peaks.DS[[4]] - min(peaks.DS[[4]]))/(max(peaks.DS[[4]]) - min(peaks.DS[[4]]))
peaks.NC[[4]] <- (peaks.NC[[4]] - min(peaks.NC[[4]]))/(max(peaks.NC[[4]]) - min(peaks.NC[[4]]))
peaks.RE[[4]] <- (peaks.RE[[4]] - min(peaks.RE[[4]]))/(max(peaks.RE[[4]]) - min(peaks.RE[[4]]))

# filter cytoband data
chr.include <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12")
cyto <- cyto[cyto[[1]] %in% chr.include, ] #exclue unwanted contigs or scaffolds
layout_nrow <- length(unique(cyto[[1]]))  #set layout height

#Get Promoter-Promoter Interacting clusters
H3K9ac.DS.PPI.bedpe <- anno.clust.H3K9ac.DS%>%
  filter(annotation1=="Promoter",annotation2=="Promoter",type==1)%>%
  dplyr::select(anchor1.chr, anchor1.str, anchor1.end,anchor2.chr, anchor2.str, anchor2.end,ipet.counts)%>%unique()

H3K9ac.NC.PPI.bedpe <- anno.clust.H3K9ac.NC%>%
  filter(annotation1=="Promoter",annotation2=="Promoter",type==1)%>%
  dplyr::select(anchor1.chr, anchor1.str, anchor1.end,anchor2.chr, anchor2.str, anchor2.end,ipet.counts)%>%unique()

H3K9ac.RE.PPI.bedpe <- anno.clust.H3K9ac.RE%>%
  filter(annotation1=="Promoter",annotation2=="Promoter",type==1)%>%
  dplyr::select(anchor1.chr, anchor1.str, anchor1.end,anchor2.chr, anchor2.str, anchor2.end,ipet.counts)%>%unique()

#plot loop map
pdf(file = "../05_General_analysis_based_on_clusters/00_Interacting_profiles/06_H3K9ac_NC.Loopgraph.pdf", width = 8, height = 10,bg=bg)
#Set ploting layout
plot_layout(nrow = layout_nrow * 3,
            ncol = 2,
            heights = unit(rep(c(1.3, 0.4, 1.3)/(layout_nrow * 3), layout_nrow), "npc"),
            widths = unit(c(0.1, 0.9), "npc"))

plot_intra_chr_interaction(interaction = clust.H3K9ac.NC%>%filter(type==1)%>%select(chrom1,start1,end1,chrom2,start2,end2,ipet.counts),
                           cyto = cyto, peaks = peaks.NC)
dev.off()

pdf(file = "../05_General_analysis_based_on_clusters/00_Interacting_profiles/07_H3K9ac_DS.Loopgraph.pdf", width = 8, height = 10,bg=bg)
#Set ploting layout
plot_layout(nrow = layout_nrow * 3,
            ncol = 2,
            heights = unit(rep(c(1.3, 0.4, 1.3)/(layout_nrow * 3), layout_nrow), "npc"),
            widths = unit(c(0.1, 0.9), "npc"))

plot_intra_chr_interaction(interaction = clust.H3K9ac.DS%>%filter(type==1)%>%select(chrom1,start1,end1,chrom2,start2,end2,ipet.counts),
                           cyto = cyto, peaks = peaks.DS)
dev.off()

pdf(file = "../05_General_analysis_based_on_clusters/00_Interacting_profiles/08_H3K9ac_RE.Loopgraph.pdf", width = 8, height = 10,bg=bg)
#Set ploting layout
plot_layout(nrow = layout_nrow * 3,
            ncol = 2,
            heights = unit(rep(c(1.3, 0.4, 1.3)/(layout_nrow * 3), layout_nrow), "npc"),
            widths = unit(c(0.1, 0.9), "npc"))

plot_intra_chr_interaction(interaction = clust.H3K9ac.RE%>%filter(type==1)%>%select(chrom1,start1,end1,chrom2,start2,end2,ipet.counts),
                           cyto = cyto, peaks = peaks.RE)
dev.off()
