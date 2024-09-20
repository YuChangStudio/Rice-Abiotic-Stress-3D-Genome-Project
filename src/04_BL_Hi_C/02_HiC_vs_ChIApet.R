library(dplyr)
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
library(Vennerable)
library(agricolae)

#color
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"

condition.col<- c("NC"=col.NC,"DS"=col.DS,"RE"=col.RE)

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
#######=========================================================================================#######
#######Assess MboI restriction fragment size.
#######=========================================================================================#######
MboI.gw <- read.delim("~/Genome_Annotation/Rice/genome_anno/MH63RS3/MH63RS3.digestMap_MboI.bed",header = F)%>%
            mutate(length=V3-V2,cat="Genome-wide")%>%select(length,cat)

MboI.H3K9ac.gw <- read.delim("../01_HPT_out/MH63RS3.MboI_RF_on_H3K9ac.bed",header = F)%>%
            mutate(length=V3-V2,cat="H3K9ac-marked regions")%>%select(length,cat)


MboI.gw$length%>%quantile()


MboI.H3K9ac.gw$length%>%quantile()

pdf(file = "./MboI_RF.length.pdf",height = 3.7,width = 3,bg=bg)
rbind(MboI.gw,MboI.H3K9ac.gw)%>%
  ggplot(aes(x=cat,y=log2(length)))+
  geom_violin()+
  geom_boxplot(width=0.4)
dev.off()



#######=========================================================================================#######
#######Load HiC loops
#######=========================================================================================#######

header.cluster <- c("chrom1","start1",	"end1",	"chrom2",	"start2",	"end2",
                    "pet.count",	"type",	"distance",	"tag.count.within.anchor.1",
                    "tag.count.within.anchor.2",	"p-value",	"p.adjust",	"-log10(p-value)",	"-log10(p.adjust)")

HIC.D.loop <- read.delim("../01_HPT_out/NEO/MH63_D_HiC_pooled_neo.cluster.FDRfiltered.txt",header = F)
HIC.N.loop <- read.delim("../01_HPT_out/NEO/MH63_N_HiC_pooled_neo.cluster.FDRfiltered.txt",header = F)
HIC.RE.loop <- read.delim("../01_HPT_out/NEO/MH63_RE_HiC_pooled_neo.cluster.FDRfiltered.txt",header = F)


colnames(HIC.D.loop) <- header.cluster
colnames(HIC.N.loop) <- header.cluster
colnames(HIC.RE.loop) <- header.cluster

pet.threshold <- 6

HIC.D.loop%<>%filter(pet.count >= pet.threshold , p.adjust < 0.05,distance > 0)%>%
  mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))%>%
  mutate(condition="DS")%>%
  select(pet.count,loop.index,condition)%>%
  filter(!(loop.index%in%ChIA.diffloop$loop.index))

HIC.N.loop%<>%filter(pet.count >= 6 , p.adjust < 0.05,distance > 0)%>%
  mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))%>%
  mutate(condition="NC")%>%
  select(pet.count,loop.index,condition)%>%
  filter(!(loop.index%in%ChIA.diffloop$loop.index))

HIC.RE.loop%<>%filter(pet.count >= pet.threshold , p.adjust < 0.05,distance > 0)%>%
  mutate(loop.index=str_c(chrom1,start1,end1,chrom2,start2,end2,sep = "_"))%>%
  mutate(condition="RE")%>%
  select(pet.count,loop.index,condition)%>%
  filter(!(loop.index%in%ChIA.diffloop$loop.index))


#######=========================================================================================#######
#######Get HiC differential loops
#######=========================================================================================#######
HIC.loops.upset <- rbind(HIC.D.loop,HIC.N.loop,HIC.RE.loop)%>%
  as_tibble()%>%
  mutate(occur=1)%>%
  tidyr::pivot_wider(
    id_cols = c("loop.index"),
    names_from = condition,
    values_from = occur,
    values_fill = list(occur = 0)
  )%>%as.data.frame()


UpSetR::upset(HIC.loops.upset,sets = c("NC","DS","RE"),
              keep.order = TRUE,
              nintersects = 7,
              order.by = "freq",
              sets.bar.color =c(col.NC,col.DS,col.RE))

#######=========================================================================================#######
#######Load ChIA-PET diffloop result
#######=========================================================================================#######
ChIA.diffloop <- read.xlsx("../../../../22_diffloop_analysis/01_diffLoops_treatments/02_diffLoop_all_list/allLoops.intersectionSet.NC.DS.RE.xlsx")


#Out put the union intersection file for Supplemental Table
ppi.loops <- read.xlsx("../../../../22_diffloop_analysis/01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")%>%pull(loop.index)
bzip23bs <- read.xlsx("../../../../25_bzip23_TFBS_MH63_ChIApet_integrative/bZIP23_bind_anchor_DS.xlsx")%>%select(V4,V5)

ChIA.diffloop%>%mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                left_join(HIC.loops.upset,by = "loop.index")%>%
                mutate(PPI=ifelse(loop.index %in% ppi.loops,"Yes","NO"))%>%
                left_join(bzip23bs,by=c("anchor1"="V4"))%>%dplyr::rename(bzip23Anchor1=V5)%>%
                left_join(bzip23bs,by=c("anchor2"="V4"))%>%dplyr::rename(bzip23Anchor2=V5)%>%
    write.xlsx(file = "Intersections_ChIApet_HiC.xlsx",keepNA=T,overwrite = T)

#######=========================================================================================#######
#######Compare BL-Hi-C and ChIA-PET loops under each condition.
#######=========================================================================================#######
ChIA.N.loop <- ChIA.diffloop%>%filter(loop.in.NC==1)%>%select(loop.index)
ChIA.N.loop%>%filter(loop.index%in%HIC.N.loop$loop.index)%>%dim() ##intersect/ChIA=19199/24235=0.7922014



ChIA.D.loop <- ChIA.diffloop%>%filter(loop.in.DS==1)%>%select(loop.index)
ChIA.D.loop%>%filter(loop.index%in%HIC.D.loop$loop.index)%>%dim() ##intersect/ChIA=9204/9577=0.9610525


ChIA.RE.loop <- ChIA.diffloop%>%filter(loop.in.RE==1)%>%select(loop.index)
ChIA.RE.loop%>%filter(loop.index%in%HIC.RE.loop$loop.index)%>%dim() ##intersect/ChIA=8626/10626=0.8117824

#plot overlapped histogram
pdf(file = "Percentage_ChIAloop_by_HIC.pdf",height = 3,width = 3,bg=bg)
data.frame(loop.count=c(24235,19199,9577,9204,10626,8626),cat=c("01CHIA_NC","02HIC_NC","03CHIA_DS","04HIC_DS","05CHIA_RE","06HIC_RE"))%>%
  ggplot(.,aes(x=cat,y=loop.count))+
  geom_bar(stat = 'identity')
dev.off()

#plot veen
pdf(file = "ChIAloop_olp_HIC.NC.veen.pdf",height = 3,width = 3,bg=bg)
Vennerable::Venn(Sets = list("ChIA-PET"=ChIA.N.loop$loop.index,"BL-Hi-C"=HIC.N.loop$loop.index))%>%plot()
dev.off()

pdf(file = "ChIAloop_olp_HIC.DS.veen.pdf",height = 3,width = 3,bg=bg)
Vennerable::Venn(Sets = list("ChIA-PET"=ChIA.D.loop$loop.index,"BL-Hi-C"=HIC.D.loop$loop.index))%>%plot()
dev.off()

pdf(file = "ChIAloop_olp_HIC.RE.veen.pdf",height = 3,width = 3,bg=bg)
Vennerable::Venn(Sets = list("ChIA-PET"=ChIA.RE.loop$loop.index,"BL-Hi-C"=HIC.RE.loop$loop.index))%>%plot()
dev.off()

#######=========================================================================================#######
#######Compare the diffloops by BL-Hi-C and ChIA-PET.
#######=========================================================================================#######
#NC->DS
ChIA.lost.NC_DS <- ChIA.diffloop%>%filter(loop.in.NC==1,loop.in.DS==0)
HIC.lost.NC_DS <- HIC.loops.upset%>%filter(NC==1,DS==0)
ChIA.lost.NC_DS%>%filter(loop.index%in%HIC.lost.NC_DS$loop.index)%>%dim() #12004/16989=0.7065748


ChIA.gain.NC_DS <- ChIA.diffloop%>%filter(loop.in.DS==1,loop.in.NC==0)
HIC.gain.NC_DS <- HIC.loops.upset%>%filter(NC==0,DS==1)
ChIA.gain.NC_DS%>%filter(loop.index%in%HIC.gain.NC_DS$loop.index)%>%dim() #2064/2331=0.8854569


ChIA.const.NC_DS <- ChIA.diffloop%>%filter(loop.in.DS==1,loop.in.NC==1)
HIC.const.NC_DS <- HIC.loops.upset%>%filter(NC==1,DS==1)
ChIA.const.NC_DS%>%filter(loop.index%in%HIC.const.NC_DS$loop.index)%>%dim() #7131/7246=0.9841292

#DS->RE
ChIA.lost.DS_RE <- ChIA.diffloop%>%filter(loop.in.DS==1,loop.in.RE==0)
HIC.lost.DS_RE <- HIC.loops.upset%>%filter(DS==1,RE==0)
ChIA.lost.DS_RE%>%filter(loop.index%in%HIC.lost.DS_RE$loop.index)%>%dim() #4061/4387=0.9256895


ChIA.gain.DS_RE <- ChIA.diffloop%>%filter(loop.in.DS==0,loop.in.RE==1)
HIC.gain.DS_RE <- HIC.loops.upset%>%filter(DS==0,RE==1)
ChIA.gain.DS_RE%>%filter(loop.index%in%HIC.gain.DS_RE$loop.index)%>%dim() #3477/5436=0.6396247

ChIA.const.DS_RE <- ChIA.diffloop%>%filter(loop.in.DS==1,loop.in.RE==1)
HIC.const.DS_RE <- HIC.loops.upset%>%filter(DS==1,RE==1)
ChIA.const.DS_RE%>%filter(loop.index%in%HIC.const.DS_RE$loop.index)%>%dim() #5128/5190=0.9880539

#plot overlapped histogram
pdf(file = "Percentage_diffloop_by_HIC.pdf",height = 4,width = 4,bg=bg)
data.frame(loop.count=c(16989,12004,2331,2064,7246,7131,4387,4061,5436,3477,5190,5128),
           cat=c("01CHIA_LOST_ND","02HIC_LOST_ND","03CHIA_GAIN_ND","04HIC_GAIN_ND","05CHIA_const_ND","06HIC_const_ND",
                 "07CHIA_LOST_DR","08HIC_LOST_DR","09CHIA_GAIN_DR","10HIC_GAIN_DR","11CHIA_const_DR","12HIC_const_DR"))%>%
  ggplot(.,aes(x=cat,y=loop.count))+
  geom_bar(stat = 'identity')+
  coord_flip()
dev.off()

#######=========================================================================================#######
#######Compare H3K9ac change levels on the contrast loops by HIC and CHIAPET
#######=========================================================================================#######
DN_quant <-read.xlsx("../../../../01_ChIPseq_MH63/04_DBA/02_spreadsheet/dba.all.MH_D_H3K9ac_vs_MH_N_H3K9ac.xlsx")%>%
            mutate(anchor.id=str_c(seqnames,start,end,sep="_"))%>%
            select(anchor.id,Conc_MH_D_H3K9ac,Conc_MH_N_H3K9ac,Fold,FDR)

DR_quant <-read.xlsx("../../../../01_ChIPseq_MH63/04_DBA/02_spreadsheet/dba.all.MH_RE_H3K9ac_vs_MH_D_H3K9ac.xlsx")%>%
            mutate(anchor.id=str_c(seqnames,start,end,sep="_"))%>%
            select(anchor.id,Conc_MH_RE_H3K9ac,Conc_MH_D_H3K9ac,Fold,FDR)


ChIAHIC.lost.NC_DS.anchor <-c(ChIA.lost.NC_DS%>%filter(loop.index%in%HIC.lost.NC_DS$loop.index)%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor1),
                           ChIA.lost.NC_DS%>%filter(loop.index%in%HIC.lost.NC_DS$loop.index)%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor2))%>%unique()

ChIAHIC.gain.NC_DS.anchor <-c(ChIA.gain.NC_DS%>%filter(loop.index%in%HIC.gain.NC_DS$loop.index)%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor1),
                           ChIA.gain.NC_DS%>%filter(loop.index%in%HIC.gain.NC_DS$loop.index)%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor2))%>%unique()

ChIAHIC.const.NC_DS.anchor <-c(ChIA.const.NC_DS%>%filter(loop.index%in%HIC.const.NC_DS$loop.index)%>%
                                 mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                        anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                 mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                 pull(anchor1),
                               ChIA.const.NC_DS%>%filter(loop.index%in%HIC.const.NC_DS$loop.index)%>%
                                 mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                        anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                 mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                 pull(anchor2))%>%unique()

ChIA.lost.NC_DS.anchor <-c(ChIA.lost.NC_DS%>%filter(!(loop.index%in%HIC.lost.NC_DS$loop.index))%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor1),
                              ChIA.lost.NC_DS%>%filter(!(loop.index%in%HIC.lost.NC_DS$loop.index))%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor2))%>%unique()

ChIA.gain.NC_DS.anchor <-c(ChIA.gain.NC_DS%>%filter(!(loop.index%in%HIC.gain.NC_DS$loop.index))%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor1),
                              ChIA.gain.NC_DS%>%filter(!(loop.index%in%HIC.gain.NC_DS$loop.index))%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor2))%>%unique()

ChIA.const.NC_DS.anchor <-c(ChIA.const.NC_DS%>%filter(!(loop.index%in%HIC.const.NC_DS$loop.index))%>%
                              mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                     anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                              mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                              pull(anchor1),
                            ChIA.const.NC_DS%>%filter(!(loop.index%in%HIC.const.NC_DS$loop.index))%>%
                              mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                     anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                              mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                              pull(anchor2))%>%unique()



ChIAHIC.lost.DS_RE.anchor <-c(ChIA.lost.DS_RE%>%filter(loop.index%in%HIC.lost.DS_RE$loop.index)%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor1),
                              ChIA.lost.DS_RE%>%filter(loop.index%in%HIC.lost.DS_RE$loop.index)%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor2))%>%unique()

ChIAHIC.gain.DS_RE.anchor <-c(ChIA.gain.DS_RE%>%filter(loop.index%in%HIC.gain.DS_RE$loop.index)%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor1),
                              ChIA.gain.DS_RE%>%filter(loop.index%in%HIC.gain.DS_RE$loop.index)%>%
                                mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                       anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                pull(anchor2))%>%unique()

ChIAHIC.const.DS_RE.anchor <-c(ChIA.const.DS_RE%>%filter(loop.index%in%HIC.const.DS_RE$loop.index)%>%
                                 mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                        anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                 mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                 pull(anchor1),
                               ChIA.const.DS_RE%>%filter(loop.index%in%HIC.const.DS_RE$loop.index)%>%
                                 mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                        anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                                 mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                                 pull(anchor2))%>%unique()

ChIA.lost.DS_RE.anchor <-c(ChIA.lost.DS_RE%>%filter(!(loop.index%in%HIC.lost.DS_RE$loop.index))%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor1),
                           ChIA.lost.DS_RE%>%filter(!(loop.index%in%HIC.lost.DS_RE$loop.index))%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor2))%>%unique()

ChIA.gain.DS_RE.anchor <-c(ChIA.gain.DS_RE%>%filter(!(loop.index%in%HIC.gain.DS_RE$loop.index))%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor1),
                           ChIA.gain.DS_RE%>%filter(!(loop.index%in%HIC.gain.DS_RE$loop.index))%>%
                             mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                    anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                             mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                             pull(anchor2))%>%unique()

ChIA.const.DS_RE.anchor <-c(ChIA.const.DS_RE%>%filter(!(loop.index%in%HIC.const.DS_RE$loop.index))%>%
                              mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                     anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                              mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                              pull(anchor1),
                            ChIA.const.DS_RE%>%filter(!(loop.index%in%HIC.const.DS_RE$loop.index))%>%
                              mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
                                     anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
                              mutate(anchor2=str_c("C",anchor2,sep=""))%>%
                              pull(anchor2))%>%unique()



##ND
#LOST
ChIAHIC.lost.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIAHIC.lost.NC_DS.anchor)%>%
                              select(Conc_MH_D_H3K9ac)%>%
                              mutate(group="ChIAHIC.lost.NC_DS",cat="lost",olp="ChIAHIC")

ChIA.lost.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIA.lost.NC_DS.anchor)%>%
                              select(Conc_MH_D_H3K9ac)%>%
                              mutate(group="ChIA.lost.NC_DS",cat="lost",olp="ChIA")

#GAIN
ChIAHIC.gain.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIAHIC.gain.NC_DS.anchor)%>%
  select(Conc_MH_D_H3K9ac)%>%
  mutate(group="ChIAHIC.gain.NC_DS",cat="gain",olp="ChIAHIC")

ChIA.gain.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIA.gain.NC_DS.anchor)%>%
  select(Conc_MH_D_H3K9ac)%>%
  mutate(group="ChIA.gain.NC_DS",cat="gain",olp="ChIA")

#CONST
ChIAHIC.const.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIAHIC.const.NC_DS.anchor)%>%
  select(Conc_MH_D_H3K9ac)%>%
  mutate(group="ChIAHIC.const.NC_DS",cat="const",olp="ChIAHIC")

ChIA.const.NC_DS.H3K9ac <- DN_quant%>%filter(anchor.id %in% ChIA.const.NC_DS.anchor)%>%
  select(Conc_MH_D_H3K9ac)%>%
  mutate(group="ChIA.const.NC_DS",cat="const",olp="ChIA")


anchor_H3K9ac_ND_plotin <- rbind(ChIAHIC.lost.NC_DS.H3K9ac,ChIA.lost.NC_DS.H3K9ac,
                                  ChIAHIC.gain.NC_DS.H3K9ac,ChIA.gain.NC_DS.H3K9ac,
                                  ChIAHIC.const.NC_DS.H3K9ac,ChIA.const.NC_DS.H3K9ac)

aov(data = anchor_H3K9ac_ND_plotin,Conc_MH_D_H3K9ac~group)%>%HSD.test(.,"group",console=TRUE)
lm(data = anchor_H3K9ac_ND_plotin,Conc_MH_D_H3K9ac~group)%>%summary()

pdf("DN_anchor_H3K9ac_D_quant.pdf",height = 3,width = 4,bg=bg)
ggplot()+
  gghalves::geom_half_violin(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIAHIC"), aes(x=cat,y=Conc_MH_D_H3K9ac,fill=olp),side = "r")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_D_H3K9ac,fill=olp), side = "l")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIAHIC"),aes(x=cat,y=Conc_MH_D_H3K9ac,fill=olp), side = "r")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_D_H3K9ac,fill=olp), side = "l")+

  gghalves::geom_half_boxplot(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIAHIC"), aes(x=cat,y=Conc_MH_D_H3K9ac),side = "r",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_D_H3K9ac), side = "l",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIAHIC"),aes(x=cat,y=Conc_MH_D_H3K9ac), side = "r",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_ND_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_D_H3K9ac), side = "l",width=0.3, outlier.shape = NA)+

  coord_flip()
dev.off()

##DR
#LOST
ChIAHIC.lost.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIAHIC.lost.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIAHIC.lost.DS_RE",cat="lost",olp="ChIAHIC")

ChIA.lost.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIA.lost.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIA.lost.DS_RE",cat="lost",olp="ChIA")

#GAIN
ChIAHIC.gain.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIAHIC.gain.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIAHIC.gain.DS_RE",cat="gain",olp="ChIAHIC")

ChIA.gain.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIA.gain.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIA.gain.DS_RE",cat="gain",olp="ChIA")

#CONST
ChIAHIC.const.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIAHIC.const.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIAHIC.const.DS_RE",cat="const",olp="ChIAHIC")

ChIA.const.DS_RE.H3K9ac <- DR_quant%>%filter(anchor.id %in% ChIA.const.DS_RE.anchor)%>%
  select(Conc_MH_RE_H3K9ac)%>%
  mutate(group="ChIA.const.DS_RE",cat="const",olp="ChIA")


anchor_H3K9ac_DR_plotin <- rbind(ChIAHIC.lost.DS_RE.H3K9ac,ChIA.lost.DS_RE.H3K9ac,
                                 ChIAHIC.gain.DS_RE.H3K9ac,ChIA.gain.DS_RE.H3K9ac,
                                 ChIAHIC.const.DS_RE.H3K9ac,ChIA.const.DS_RE.H3K9ac)

aov(data = anchor_H3K9ac_DR_plotin,Conc_MH_RE_H3K9ac~group)%>%HSD.test(.,"group",console=TRUE)


pdf("DR_anchor_H3K9ac_RE_quant.pdf",height = 3,width = 4,bg=bg)
ggplot()+
  gghalves::geom_half_violin(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIAHIC"), aes(x=cat,y=Conc_MH_RE_H3K9ac,fill=olp),side = "r")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_RE_H3K9ac,fill=olp), side = "l")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIAHIC"),aes(x=cat,y=Conc_MH_RE_H3K9ac,fill=olp), side = "r")+
  gghalves::geom_half_violin(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_RE_H3K9ac,fill=olp), side = "l")+

  gghalves::geom_half_boxplot(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIAHIC"), aes(x=cat,y=Conc_MH_RE_H3K9ac),side = "r",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_RE_H3K9ac), side = "l",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIAHIC"),aes(x=cat,y=Conc_MH_RE_H3K9ac), side = "r",width=0.3, outlier.shape = NA)+
  gghalves::geom_half_boxplot(data = anchor_H3K9ac_DR_plotin%>%filter(olp=="ChIA"),aes(x=cat,y=Conc_MH_RE_H3K9ac), side = "l",width=0.3, outlier.shape = NA)+

  coord_flip()
dev.off()



######
######network and CIA analysis
library(igraph)
library(tidygraph)

#===
ppi.anchor <- read.xlsx("../../../../20_ChIAPET_profile_analysis/03_PPI_pairs/02_All_PPI_gene_anno/PPI_Gene.anno.allConditions.xlsx")%>%
              select(anchor.id)


HIC.D.tbl <- HIC.D.loop%>%
  mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
         anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
  mutate(anchor2=str_c("C",anchor2,sep=""))%>%
  select(anchor1,anchor2,pet.count,loop.index)%>%
  filter(anchor1%in%ppi.anchor$anchor.id,anchor2%in%ppi.anchor$anchor.id)%>%
  as_tbl_graph(PPI,directed = F)%>%
  rename(anchor.id=name)%>%
  mutate(degree=centrality_degree())%>%
  mutate(ego.group=group_components())%>%
  mutate(RW.group =group_walktrap(steps = 4,weights = pet.count))

degree.total.D <- degree(HIC.D.tbl)%>%sum()


node.info.D <- HIC.D.tbl%>%
  as.data.frame%>%
  rowwise()%>%
  mutate(anchor.range=(as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[3])-
                         as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[2])))%>%
  mutate(degree.dens = (degree * 1000 * 10000)/(degree.total.D * anchor.range))


#===
HIC.N.tbl <- HIC.N.loop%>%
  mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
         anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
  mutate(anchor2=str_c("C",anchor2,sep=""))%>%
  select(anchor1,anchor2,pet.count,loop.index)%>%
  filter(anchor1%in%ppi.anchor$anchor.id,anchor2%in%ppi.anchor$anchor.id)%>%
  as_tbl_graph(PPI,directed = F)%>%
  rename(anchor.id=name)%>%
  mutate(degree=centrality_degree())%>%
  mutate(ego.group=group_components())%>%
  mutate(RW.group =group_walktrap(steps = 4,weights = pet.count))

degree.total.N <- degree(HIC.N.tbl)%>%sum()


node.info.N <- HIC.N.tbl%>%
  as.data.frame%>%
  rowwise()%>%
  mutate(anchor.range=(as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[3])-
                         as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[2])))%>%
  mutate(degree.dens = (degree * 1000 * 10000)/(degree.total.N * anchor.range))%>%
  select(anchor.id,degree.dens)%>%
  rename(CIA.N=degree.dens)

#===
HIC.RE.tbl <- HIC.RE.loop%>%
  mutate(anchor1=str_split(loop.index,pattern = "_C",simplify = T)[,1],
         anchor2=str_split(loop.index,pattern = "_C",simplify = T)[,2])%>%
  mutate(anchor2=str_c("C",anchor2,sep=""))%>%
  select(anchor1,anchor2,pet.count,loop.index)%>%
  filter(anchor1%in%ppi.anchor$anchor.id,anchor2%in%ppi.anchor$anchor.id)%>%
  as_tbl_graph(PPI,directed = F)%>%
  rename(anchor.id=name)%>%
  mutate(degree=centrality_degree())%>%
  mutate(ego.group=group_components())%>%
  mutate(RW.group =group_walktrap(steps = 4,weights = pet.count))

degree.total.RE <- degree(HIC.RE.tbl)%>%sum()


node.info.RE <- HIC.RE.tbl%>%
  as.data.frame%>%
  rowwise()%>%
  mutate(anchor.range=(as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[3])-
                         as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[2])))%>%
  mutate(degree.dens = (degree * 1000 * 10000)/(degree.total.RE * anchor.range))%>%
  select(anchor.id,degree.dens)%>%
  rename(CIA.RE=degree.dens)
#===
#Load ChIA-PET SPRs
DS.specific.SPR <- read.xlsx("../../../../23_PPI_network_and_enhancer_analysis/01_analysis_on_MH63WT/03_Super_Anchor/03_Dominate_superAnchor_list/DS_dominate_superAnchor.xlsx")

df <- node.info.D%>%filter(anchor.id %in% DS.specific.SPR$anchor.id)%>%
  select(anchor.id,degree.dens)%>%
  rename(CIA.D=degree.dens)%>%
  left_join(.,node.info.N)%>%
  left_join(.,node.info.RE)%>%
  as.data.frame()

df[is.na(df)] <-0.01

df<-df%>%mutate(D_vs_N=CIA.D/CIA.N,D_vs_RE=CIA.D/CIA.RE)%>%
          mutate(N_SIG = ifelse(D_vs_N > 1,"Yes","NO"))%>%
          mutate(RE_SIG = ifelse(D_vs_RE > 1,"Yes","NO"))
pdf("SPR.heatmap.pdf",height = 5,width = 5,bg=bg)
pheatmap::pheatmap(df%>%select(-anchor.id,-D_vs_N,-D_vs_RE,-N_SIG,-RE_SIG)%>%as.matrix(rownames.force = TRUE),
                   col=viridis::inferno(n = 255,alpha = 1,direction = -1),
                   annotation_row=df%>%select(N_SIG,RE_SIG),
                   clustering_method = "complete",
                   annotation_colors  = list(
                     N_SIG=c(Yes=col.NC,NO="white"),
                     RE_SIG=c(Yes=col.RE,NO="white"))
)
dev.off()

df%>%filter(D_vs_N < 1 ,D_vs_RE < 1)%>%dim()
df%>%filter((D_vs_N < 1 | D_vs_RE < 1))%>%dim()
df%>%filter(D_vs_N > 1 ,D_vs_RE > 1)%>%dim()

write.xlsx(x = df,file = "HIC_SPR.xlsx")
