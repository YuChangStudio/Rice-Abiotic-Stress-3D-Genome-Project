#########===============================================================================================#########
#########In this section, we analyze the occupation of ChIA-PET loop anchors regions by bZIP23 binding
#########===============================================================================================#########
library(openxlsx)
library(ggplot2)
library(ral)
library(UpSetR)
library(multcomp)
library(stringr)
library(dplyr)
library(bedtoolsr)

options(bedtools.path = "/Users/juliuschiong/opt/anaconda3/bin/") #mac
#options(bedtools.path = "wsl /mnt/d/linux_utility/anaconda3/bin") #win; result file should be read back manually from temp files under win.

#color preset
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"
col.bzip23 <- ralcolors["RAL6027"]%>%as.vector #"#81C0BB"

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
      legend.title = element_text(colour='black',size = 8,face = "bold"),
      legend.background = element_rect(fill = "white",colour = "black"),
      legend.key = element_blank())%>%theme_set()

#Load anchor list
total.anchors <- read.delim("../../01_ChIPseq_MH63/01_Peak_set/04_Unified/MH_H3K9ac_unify.IDR.narrowPeak",header = F)%>%
  mutate(anchor.id=str_c(V1,V2,V3,sep = "_"))%>%
  rowwise()%>%
  mutate(start=V2-300,end=V3+300)%>% #span regions by 300 bp (two nucleosomes' range)
  rename(chr=V1)%>%
  select(chr,start,end,anchor.id)%>%
  as.data.frame()


#Load bZIP23 DS peak regions and intensities
bzip23.IDR.DS <- read.delim("../../04_ChIPseq_anti_bZIP23/05_peaks/04_IDR/DS_bzip23.uniq.bed",header = F)

#########===============================================================================================#########
#########Find intersections between total anchors and bZIP23 TFBS
#########===============================================================================================#########
#Use the IDR peak regions under DS for overlap analysis.

bt.intersect(a = total.anchors,b=bzip23.IDR.DS,c=T)%>%
  filter(V5>0)%>%
  write.xlsx("../bZIP23_bind_anchor_DS.xlsx")


#########===============================================================================================#########
#########Find intersections DS dominate anchors and bZIP23 TFBS
#########===============================================================================================#########
bZIP23_bind_anchor.DS <- read.xlsx("../bZIP23_bind_anchor_DS.xlsx")%>%pull(V4)

#Find how many DS-dominate loops are with at least one anchor bound by bzip23
DS.dominate <- read.xlsx("../../22_diffloop_analysis/01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")%>%
                filter(loop.in.DS==1,loop.in.NC==0,loop.in.RE==0)
DS.dominate%>%filter(anchor1.id %in% bZIP23_bind_anchor.DS | anchor2.id %in% bZIP23_bind_anchor.DS)%>%dim() #1080 vs 1432; 0.7541899

#Find how many DS-dominate loops that lost in osbzip23 mutant are with at least one anchor bound by bzip23
DS.dominate.lost <- read.xlsx("../../22_diffloop_analysis/02_diffLoops_MH63_vs_bzip23/04_diffPPI_list/DS_vs_bzip23.PPI.intersectionSet.xlsx")%>%
                    filter(loop.in.DS==1,loop.in.bzip23==0,loop.index%in%DS.dominate$loop.index) #1179
DS.dominate.lost%>%filter(anchor1.id %in% bZIP23_bind_anchor.DS | anchor2.id %in% bZIP23_bind_anchor.DS)%>%dim() #881 vs 1179; 0.7472434


#Find how many DS-PPIs are with at least one anchor bound by bzip23
DS.PPI <- read.xlsx("../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/D.PPI.uniqGene.xlsx")

DS.PPI%>%filter(anchor1.id %in% bZIP23_bind_anchor.DS | anchor2.id %in% bZIP23_bind_anchor.DS)%>%dim() #5223 vs 7059; 0.7399065

pdf("../proportion_bzip23_binding_PPI_DS.pdf",height = 3,width = 8,bg=bg)
data.frame(PPI.category=c("DS PPI","DS dominate PPI","DS dominate PPI lost in osbzip23","DS PPI","DS dominate PPI","DS dominate PPI lost in osbzip23"),
           count=c(5223,1080,881,7059,1432,1179),
           group=c("bzip23_bind","bzip23_bind","bzip23_bind","total","total","total"))%>%
ggplot(aes(y=PPI.category,x=count,group=group))+
  geom_bar(position = "dodge",stat = "identity",aes(fill=group))
dev.off()



