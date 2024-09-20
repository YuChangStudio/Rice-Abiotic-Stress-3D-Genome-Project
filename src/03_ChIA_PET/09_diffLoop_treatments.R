####=========================================================================================================######
####Load overall configs; in this section, we only analyze NC, DS, and RE samples.
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
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"

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

NC.loop.all <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/N.all.loop.bedpe",header = F,col.names=loop.header)%>%
                mutate(distance=ifelse(chr.1==chr.2,round((end.2+start.2)/2-(end.1+start.1)/2),NA))%>%
                select(loop.index,distance)%>%mutate(condition="NC")

RE.loop.all <- read.delim("../../../20_ChIAPET_profile_analysis/01_loop_BEDPE/RE.all.loop.bedpe",header = F,col.names=loop.header)%>%
                mutate(distance=ifelse(chr.1==chr.2,round((end.2+start.2)/2-(end.1+start.1)/2),NA))%>%
                select(loop.index,distance)%>%mutate(condition="RE")


#PPI loops will be highlighted in the diff loops, and the TPM will be used for subsequent evaluation:
#Load only the unique gene files, the file contains the average TPM of genes in the anchor, which makes the comparision easier.
#Exclude the ipet.counts values in each file, since the value can be different for each file, and may lead to duplications in merge.

PPI.NC <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/N.PPI.uniqGene.xlsx")%>%
          mutate(condition="NC")%>%
          select(-ipet.counts)

PPI.DS <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/D.PPI.uniqGene.xlsx")%>%
          mutate(condition="DS")%>%
          select(-ipet.counts)

PPI.RE <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/RE.PPI.uniqGene.xlsx")%>%
          mutate(condition="RE")%>%
          select(-ipet.counts)

###============================================================================================###
##Analyze the overlaping/diff interaction of all types of loops under NC DS RE
###============================================================================================###
##Plot the upset plots for all loops and PPI loops separately and merge them in AI
NC.DS.RE.all.loops.upset <- rbind(DS.loop.all,NC.loop.all,RE.loop.all)%>%
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

pdf("../01_diifLoop_all_plot/All_loops.upset.pdf",height = 2.8,width = 3.5,bg=bg)
upset(NC.DS.RE.all.loops.upset,sets = c("NC","DS","RE"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      sets.bar.color =c(col.NC,col.DS,col.RE))
dev.off()

#Export the loop list for each set
rbind(NC.DS.RE.all.loops.upset%>%filter(NC==1,DS==0,RE==0)%>%mutate(set="set1"),
      NC.DS.RE.all.loops.upset%>%filter(NC==1,DS==1,RE==1)%>%mutate(set="set2"),
      NC.DS.RE.all.loops.upset%>%filter(NC==0,DS==0,RE==1)%>%mutate(set="set3"),
      NC.DS.RE.all.loops.upset%>%filter(NC==1,DS==0,RE==1)%>%mutate(set="set4"),
      NC.DS.RE.all.loops.upset%>%filter(NC==1,DS==1,RE==0)%>%mutate(set="set5"),
      NC.DS.RE.all.loops.upset%>%filter(NC==0,DS==1,RE==0)%>%mutate(set="set6"),
      NC.DS.RE.all.loops.upset%>%filter(NC==0,DS==1,RE==1)%>%mutate(set="set7"))%>%
  select(loop.index,NC,DS,RE,set)%>%
  rename(loop.in.NC=NC,loop.in.DS=DS,loop.in.RE=RE)%>%
  write.xlsx("../02_diffLoop_all_list/allLoops.intersectionSet.NC.DS.RE.xlsx")

##============================================================================##
##Analyze the overlaping/diff interaction of PPI loops under NC DS RE
##============================================================================##
#Note the expression level/FC of a PPI is expressed as the average TPM/FC of the anchor pair.
id_cols = c("loop.index", "anchor1.id","anchor2.id",
            "log2p1.DS.TPM.loop","log2p1.NC.TPM.loop","log2p1.RE.TPM.loop","log2p1.bzip23.TPM.loop",
            "log2FC_DS_vs_NC.loop","log2FC_RE_vs_DS.loop","log2FC_RE_vs_NC.loop")

NC.DS.RE.PPI.upset <- rbind(PPI.DS,PPI.NC,PPI.RE)%>%
  mutate(log2p1.DS.TPM.loop=log2((DS.TPM.anchor.1+DS.TPM.anchor.2)/2+1),
         log2p1.NC.TPM.loop=log2((NC.TPM.anchor.1+NC.TPM.anchor.2)/2+1),
         log2p1.RE.TPM.loop=log2((RE.TPM.anchor.1+RE.TPM.anchor.2)/2+1),
         log2p1.bzip23.TPM.loop=log2((bzip23.TPM.anchor.1+bzip23.TPM.anchor.2)/2+1),
         log2FC_DS_vs_NC.loop=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2,
         log2FC_RE_vs_DS.loop=(log2FC_RE_vs_DS.anchor.1+log2FC_RE_vs_DS.anchor.2)/2,
         log2FC_RE_vs_NC.loop=(log2FC_RE_vs_NC.anchor.1+log2FC_RE_vs_NC.anchor.2)/2)%>%
  as_tibble()%>%
  mutate(occur=1)%>%
  pivot_wider(
    id_cols = id_cols,
    names_from = condition,
    values_from = occur,
    values_fill = list(occur = 0)
  )%>%as.data.frame()


pdf("../03_diffPPI_plot/PPI_loops.upset.pdf",height = 2.8,width = 3.5,bg=bg)
upset(NC.DS.RE.PPI.upset,sets = c("NC","DS","RE"),
      keep.order = TRUE,
      nintersects = 7,
      order.by = "freq",
      sets.bar.color =c(col.NC,col.DS,col.RE))
dev.off()



#Plot the TPM / FC values separately by ggplot2. The boxplot.summary in upsetR is full of bugs.
#The set number is in accordence to the corresponding upset plot.
NC.DS.RE.PPI.upset.data <- rbind(NC.DS.RE.PPI.upset%>%filter(NC==1,DS==0,RE==0)%>%mutate(set="set1"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==1,DS==1,RE==1)%>%mutate(set="set2"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==0,DS==0,RE==1)%>%mutate(set="set3"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==1,DS==0,RE==1)%>%mutate(set="set4"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==1,DS==1,RE==0)%>%mutate(set="set5"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==0,DS==1,RE==0)%>%mutate(set="set6"),
                                 NC.DS.RE.PPI.upset%>%filter(NC==0,DS==1,RE==1)%>%mutate(set="set7"))%>%
  select(set,log2p1.DS.TPM.loop,log2p1.NC.TPM.loop,log2p1.RE.TPM.loop,log2p1.bzip23.TPM.loop,
         log2FC_DS_vs_NC.loop,log2FC_RE_vs_DS.loop,log2FC_RE_vs_NC.loop)%>%
  reshape2::melt(id.var = "set")

NC.DS.RE.PPI.upset.data[is.na(NC.DS.RE.PPI.upset.data)]<-0

pdf("../03_diffPPI_plot/PPI_loops.TPM.FC.boxplot.pdf",height = 12,width = 6,bg=bg)
ggplot(NC.DS.RE.PPI.upset.data,aes(x=set,y=value))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(facets = .~variable,ncol = 2,scales = "free")+
  ylab("log2(TPM+1)/log2FC")
dev.off()

#Anova analysis for each category:
detach("package:dplyr")
library(multcomp)
library(dplyr)

NC.DS.RE.PPI.upset.data$variable%>%unique()%>%
lapply(., function(x){

  df <-NC.DS.RE.PPI.upset.data%>%
    filter(variable==paste0(x))%>%
    mutate(set=as.factor(set))%>%
    lm(value ~ set,data = .)%>%
    glht(.,linfct = mcp(set = "Tukey"))%>%cld(decreasing=T)

  df$mcletters$Letters%>%as.data.frame()%>%rename(hsd=".")%>%mutate(set=rownames(.))%>%
    assign(paste0("NC.DS.RE.PPI.upset.",x,".hsd"),.,envir = .GlobalEnv)
})


#Export the PPI loop information of each intersection set.

rbind(NC.DS.RE.PPI.upset%>%filter(NC==1,DS==0,RE==0)%>%mutate(set="set1"),
      NC.DS.RE.PPI.upset%>%filter(NC==1,DS==1,RE==1)%>%mutate(set="set2"),
      NC.DS.RE.PPI.upset%>%filter(NC==0,DS==0,RE==1)%>%mutate(set="set3"),
      NC.DS.RE.PPI.upset%>%filter(NC==1,DS==0,RE==1)%>%mutate(set="set4"),
      NC.DS.RE.PPI.upset%>%filter(NC==1,DS==1,RE==0)%>%mutate(set="set5"),
      NC.DS.RE.PPI.upset%>%filter(NC==0,DS==1,RE==0)%>%mutate(set="set6"),
      NC.DS.RE.PPI.upset%>%filter(NC==0,DS==1,RE==1)%>%mutate(set="set7"))%>%
  select(loop.index,NC,DS,RE,set)%>%
  rename(loop.in.NC=NC,loop.in.DS=DS,loop.in.RE=RE)%>%
  left_join(.,rbind(PPI.DS,PPI.NC,PPI.RE)%>%select(-condition)%>%unique(),by="loop.index")%>%
  write.xlsx("../04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")