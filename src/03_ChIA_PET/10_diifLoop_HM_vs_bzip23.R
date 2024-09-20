####=========================================================================================================######
####Load overall configs; in this section, we only analyze MH63-WT and osbzip23 samples.
####=========================================================================================================######
library(stringr)
library(openxlsx)
library(ggplot2)
library(ral)
library(tidyr)
library(UpSetR)
library(dplyr)
library(magrittr)

#color preset
bg="transparent"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.bzip23 <- ralcolors["RAL6027"]%>%as.vector #"#81C0BB"

#ggplot2 theme
theme(panel.background = element_rect(fill = "transparent",colour='black',size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
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

#Load loop.
#All categories of loops are load for first-step comparison.

loop.header <- c("chr.1","start.1","end.1","chr.2","start.2","end.2","loop.index","ipet.count")

DS.loop.all <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/D.all.loop.bedpe",header = F,col.names=loop.header)%>%
                mutate(distance=ifelse(chr.1==chr.2,round((end.2+start.2)/2-(end.1+start.1)/2),NA))%>%
                select(loop.index,distance)%>%mutate(condition="DS")


bzip23.loop.all <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/bzip23.all.loop.bedpe",header = F,col.names=loop.header)%>%
                mutate(distance=ifelse(chr.1==chr.2,round((end.2+start.2)/2-(end.1+start.1)/2),NA))%>%
                select(loop.index,distance)%>%mutate(condition="bzip23")

#PPI loops will be highlighted in the diff loops, and the TPM will be used for subsequent evaluation:
#Load only the unique gene files, the file contains the average TPM of genes in the anchor, which makes the comparision easier.

PPI.DS <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/D.PPI.uniqGene.xlsx")%>%mutate(condition="DS")
PPI.bzip23 <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/bzip23.PPI.uniqGene.xlsx")%>%mutate(condition="bzip23")

###============================================================================================###
##Analyze the overlaping/diff interaction of all types of loops
###============================================================================================###
##Plot the upset plots for all loops and PPI loops separately and merge them in AI
DS.bzip23.all.loops.upset <- rbind(DS.loop.all,bzip23.loop.all)%>%
                            as_tibble()%>%
                            mutate(occur=1)%>%
                            mutate(log2p1.distance=log2(distance+1))%>%
                            select(-distance)%>%
                            pivot_wider(
                              id_cols = c("loop.index","log2p1.distance"),
                              names_from = condition,
                              values_from = occur,
                              values_fill = list(occur = 0)
                            )%>%as.data.frame()

pdf("../01_diifLoop_all_plot/DS_vs_bzip23_All_loops.upset.pdf",height = 2.8,width = 3.5,bg=bg)
upset(DS.bzip23.all.loops.upset,sets = c("DS","bzip23"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      sets.bar.color =c(col.DS,col.bzip23))
dev.off()

#Export the loop list for each set
rbind(DS.bzip23.all.loops.upset%>%filter(DS==1,bzip23==0)%>%mutate(set="set1"),
      DS.bzip23.all.loops.upset%>%filter(DS==1,bzip23==1)%>%mutate(set="set2"),
      DS.bzip23.all.loops.upset%>%filter(DS==0,bzip23==1)%>%mutate(set="set3"))%>%
  select(loop.index,DS,bzip23,set)%>%
  rename(loop.in.DS=DS,loop.in.bzip23=bzip23)%>%
  write.xlsx("../02_diffLoop_all_list/DS_vs_bzip23.allLoops.intersectionSet.xlsx")

##============================================================================##
##Analyze the overlaping/diff interaction of PPI loops
##============================================================================##
#Note the expression level/FC of a PPI is expressed as the average TPM/FC of the anchor pair.
id_cols = c("loop.index", "anchor1.id","anchor2.id",
            "log2p1.DS.TPM.loop","log2p1.NC.TPM.loop","log2p1.RE.TPM.loop","log2p1.bzip23.TPM.loop",
            "log2FC_DS_vs_NC.loop","log2FC_RE_vs_DS.loop","log2FC_RE_vs_NC.loop","log2FC_DS_vs_bzip23.loop")

DS.bzip23.PPI.upset <- rbind(PPI.DS,PPI.bzip23)%>%
  mutate(log2p1.DS.TPM.loop=log2((DS.TPM.anchor.1+DS.TPM.anchor.2)/2+1),
         log2p1.NC.TPM.loop=log2((NC.TPM.anchor.1+NC.TPM.anchor.2)/2+1),
         log2p1.RE.TPM.loop=log2((RE.TPM.anchor.1+RE.TPM.anchor.2)/2+1),
         log2p1.bzip23.TPM.loop=log2((bzip23.TPM.anchor.1+bzip23.TPM.anchor.2)/2+1),
         log2FC_DS_vs_NC.loop=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2,
         log2FC_RE_vs_DS.loop=(log2FC_RE_vs_DS.anchor.1+log2FC_RE_vs_DS.anchor.2)/2,
         log2FC_RE_vs_NC.loop=(log2FC_RE_vs_NC.anchor.1+log2FC_RE_vs_NC.anchor.2)/2,
         log2FC_DS_vs_bzip23.loop=(log2FC_DS_vs_bzip23.anchor.1+log2FC_DS_vs_bzip23.anchor.2)/2)%>%
  as_tibble()%>%
  mutate(occur=1)%>%
  pivot_wider(
    id_cols = id_cols,
    names_from = condition,
    values_from = occur,
    values_fill = list(occur = 0)
  )%>%as.data.frame()


pdf("../03_diffPPI_plot/DS_vs_bzip23_PPI_loops.upset.pdf",height = 2.8,width = 3.5,bg=bg)
upset(DS.bzip23.PPI.upset,sets = c("DS","bzip23"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      sets.bar.color =c(col.DS,col.bzip23))
dev.off()

#Plot the TPM / FC values separately by ggplot2. The boxplot.summary in upsetR is full of bugs.
#The set number is in accordence to the corresponding upset plot.
DS.bzip23.PPI.upset.data <- rbind(DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==0)%>%mutate(set="set1"),
      DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==1)%>%mutate(set="set2"),
      DS.bzip23.PPI.upset%>%filter(DS==0,bzip23==1)%>%mutate(set="set3"))%>%
  select(set,log2p1.DS.TPM.loop,log2p1.NC.TPM.loop,log2p1.RE.TPM.loop,log2p1.bzip23.TPM.loop,
         log2FC_DS_vs_NC.loop,log2FC_RE_vs_DS.loop,log2FC_RE_vs_NC.loop,log2FC_DS_vs_bzip23.loop)%>%
  reshape2::melt(id.var = "set")

DS.bzip23.PPI.upset.data[is.na(DS.bzip23.PPI.upset.data )]<-0

pdf("../03_diffPPI_plot/DS_vs_bzip23.PPI_loops.TPM.FC.boxplot.pdf",height = 12,width = 6,bg=bg)
ggplot(DS.bzip23.PPI.upset.data,aes(x=set,y=value))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(facets = .~variable,ncol = 2,scales = "free")+
  ylab("log2(TPM+1)/log2FC")
dev.off()

#Anova analysis for each category:
detach("package:dplyr")
library(multcomp)
library(dplyr)

DS.bzip23.PPI.upset.data$variable%>%unique()%>%
lapply(., function(x){

  df <-DS.bzip23.PPI.upset.data%>%
    filter(variable==paste0(x))%>%
    mutate(set=as.factor(set))%>%
    lm(value ~ set,data = .)%>%
    glht(.,linfct = mcp(set = "Tukey"))%>%cld(decreasing=T)

  df$mcletters$Letters%>%as.data.frame()%>%rename(hsd=".")%>%mutate(set=rownames(.))%>%
    assign(paste0("DS.bzip23.PPI.upset.",x,".hsd"),.,envir = .GlobalEnv)
})


#Export the PPI loop information of each intersection set.
rbind(DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==0)%>%mutate(set="set1"),
      DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==1)%>%mutate(set="set2"),
      DS.bzip23.PPI.upset%>%filter(DS==0,bzip23==1)%>%mutate(set="set3"))%>%
  select(loop.index,DS,bzip23,set)%>%
  rename(loop.in.DS=DS,loop.in.bzip23=bzip23)%>%
  left_join(.,rbind(PPI.DS,PPI.bzip23)%>%select(-condition,-ipet.counts)%>%unique(),by="loop.index")%>%
  write.xlsx("../04_diffPPI_list/DS_vs_bzip23.PPI.intersectionSet.xlsx")

##See how many DS-dominate loops are in each set and label them manually in AI.
DS.dominatePPI <- read.xlsx("../../01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")%>%
                  filter(loop.in.NC==0, loop.in.DS==1, loop.in.RE==0)%>%pull(loop.index)


DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==0)%>%filter(loop.index %in% DS.dominatePPI)%>%dim() #1179
DS.bzip23.PPI.upset%>%filter(DS==1,bzip23==1)%>%filter(loop.index %in% DS.dominatePPI)%>%dim() #253
DS.bzip23.PPI.upset%>%filter(DS==0,bzip23==1)%>%filter(loop.index %in% DS.dominatePPI)%>%dim() #0
