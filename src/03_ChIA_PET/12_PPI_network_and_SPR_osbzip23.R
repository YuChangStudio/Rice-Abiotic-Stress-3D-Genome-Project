#########===============================================================================================#########
#########In this section, we analyze putative enhancer/super-promoter associating anchors by network analysis
#########===============================================================================================#########

#########===============================================================================================#########
#########Load configs
#########===============================================================================================#########
library(igraph)
library(tidygraph)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggraph)
library(ggsignif)
library(openxlsx)
library(patchwork)
library(graphlayouts)
library(ral)

#color preset
bg="transparent"
col.NC <- ralcolors["RAL6021"]%>%as.vector #"#86A47C"
col.DS <- ralcolors["RAL8002"]%>%as.vector #"#7B5141"
col.RE <- ralcolors["RAL1002"]%>%as.vector #"#D2B773"
col.bzip23 <- ralcolors["RAL6027"]%>%as.vector #"#81C0BB"

col.conditions <- c("bzip23"=col.bzip23,"DS"=col.DS)

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

#MH63 gene annotations
mh63.encode <- read.xlsx("~/genome_anno/MH63RS3/MH63RS3_encode_Sep2022.xlsx")

#########===============================================================================================#########
#########Load data for analysis on NC, DS, RE PPIs.
#########===============================================================================================#########
##We only analyze the situation of DS vs bzip23 at this time.


#Load anchor info
PPI.anchor.info <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/03_PPI_anchor_annoation/PPI_anchor_unique_genes.annoation.xlsx")%>%
                    left_join(.,mh63.encode,by=c("represent.gene"="MH63RS3_ID"))%>%
                    mutate(node.name=ifelse(name != "" , name, represent.gene))

#give 0 to the NA quantification of low expression genes
PPI.anchor.info[is.na(PPI.anchor.info)] <- 0

##Load results on condition-dominate PPIs as edge info.
PPI.condition.domination <- read.xlsx("../../../22_diffloop_analysis/01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")%>%
  mutate(loop.domination=case_when(
    (loop.in.DS ==1 & loop.in.NC==0 & loop.in.RE==0) ~ "DS_dominate",
    (loop.in.DS ==0 & loop.in.NC==1 & loop.in.RE==0) ~ "NC_dominate",
    (loop.in.DS ==0 & loop.in.NC==0 & loop.in.RE==1) ~ "RE_dominate",
    (loop.in.DS ==1 & loop.in.NC==1 & loop.in.RE==1) ~ "Conserved",
    TRUE ~ "other"
  ))%>%
  select(-anchor1.id,-anchor2.id)

PPI.condition.domination%<>%select(loop.index,loop.domination)


PPI.DS_vs_bzip23.dominate <- read.xlsx("../../../22_diffloop_analysis/02_diffLoops_MH63_vs_bzip23/04_diffPPI_list/DS_vs_bzip23.PPI.intersectionSet.xlsx")%>%
  mutate(loop.DS_vs_bzip23=case_when(
    set=="set1" ~ "lost_in_bzip23",
    set=="set2" ~ "dispensable_to_bzip23",
    set=="set3" ~ "denovo_in_bzip23"))

PPI.DS_vs_bzip23.dominate %<>% select(loop.index,loop.in.DS,loop.in.bzip23,loop.DS_vs_bzip23)

PPI.loop.class <- PPI.DS_vs_bzip23.dominate%>%left_join(PPI.condition.domination)

#########===============================================================================================#########
#########Generate tidygraph network object and make basic analysis
#########===============================================================================================#########
##Generate tbl network object:
#degree.dens=(degree*1000)/(total degree in the FG group * anchor length in bp)
#We expect that for lager group, the super anchor should have larger degree
lapply(c("D","bzip23"), function(x){

  file.prefix <- "../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/"

  #get raw edge
  #exclude anchor annoation as they will be lost during tbl generation.
  PPI <- read.xlsx(paste0(file.prefix,x,".PPI.uniqGene.xlsx"))%>%
    select(anchor1.id,anchor2.id,loop.index,ipet.counts,type)%>%#get edge info
    left_join(.,PPI.loop.class,by="loop.index")

  assign(paste0(x,".PPIloops"),PPI,envir = .GlobalEnv)

  #generate network tbl object
  tbl <- as_tbl_graph(PPI,directed = F)%>%     #anchors are used as nodes.
    rename(anchor.id=name)%>%
    mutate(degree=centrality_degree())%>% #Calculate the degree of nodes.
    mutate(ego.group=group_components())%>% #Label connected sub-networks. NO.1 is always the largest one.
    mutate(RW.group =group_walktrap(steps = 4,weights = ipet.counts), #Identify PPI communities by random walk.
           #PL.group =group_label_prop(weights = ipet.counts), #Identify PPI communities by propagating labels.
           #FG.group = group_fast_greedy(),  #Identify PPI communities by fast-greedy.
    )

  degree.total <- degree(tbl)%>%sum()

  #annotate nodes.
  node.info <- tbl%>%
    as.data.frame%>%
    left_join(PPI.anchor.info)%>%
    group_by(RW.group)%>%
    mutate(nodes.in.RWgroup=length(unique(anchor.id)), #get the total number of nodes and degrees within each fast-greedy group.
           degree.in.RWgroup=sum(degree), #summed degree within each RW groups
           avgTPM.DS.RWgroup=mean(DS.TPM.anchor), #The the average/median TPM of the fast-greedy group for subsequent filtering.
           medianTPM.DS.RWgroup=median(DS.TPM.anchor),
           avgTPM.NC.RWgroup=mean(NC.TPM.anchor),
           medianTPM.NC.RWgroup=median(NC.TPM.anchor),
           avgTPM.RE.RWgroup=mean(RE.TPM.anchor),
           medianTPM.RE.RWgroup=median(RE.TPM.anchor),
           avgTPM.bzip23.RWgroup=mean(bzip23.TPM.anchor),
           medianTPM.bzip23.RWgroup=median(bzip23.TPM.anchor))%>%
    ungroup()%>%

    group_by(ego.group)%>%
    mutate(nodes.in.ego.group=length(unique(anchor.id)), #get the total number of nodes and degrees within each ego group.
           degree.in.ego.group=sum(degree))%>%
    ungroup()%>%
    rowwise()%>%
    mutate(anchor.range=(as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[3])-
                           as.numeric(str_split(anchor.id,pattern = "_",simplify = T)[2])))%>%
    mutate(degree.dens = (degree * 1000 * 10000)/(degree.total * anchor.range)) ##Get degree.dens


  #get annotated tbl object
  tbl%>%
    left_join(.,node.info)%>%
    assign(paste0(x,".tbl"),.,envir = .GlobalEnv)
})

##RW group number under each condition


D.tbl%>%select(RW.group)%>%as.data.frame()%>%max()#1427

bzip23.tbl%>%select(RW.group)%>%as.data.frame()%>%max()#868

##Compare the RW group size (by degree):
degree.in.RWgroup.data <- rbind(
  D.tbl%>%select(degree.in.RWgroup)%>%as.data.frame()%>%mutate(condition="DS"),
  bzip23.tbl%>%select(degree.in.RWgroup)%>%as.data.frame()%>%mutate(condition="bzip23")
  #bzip23.tbl%>%select(degree.in.RWgroup)%>%as.data.frame()%>%mutate(condition="bzip23")
)

degree.in.RWgroup.data$degree.in.RWgroup <- as.numeric(degree.in.RWgroup.data$degree.in.RWgroup)
degree.in.RWgroup.data$condition <- as.factor(degree.in.RWgroup.data$condition)

##Compare the RW group size (by node number):
nodes.in.RWgroup.data <- rbind(
  D.tbl%>%select(nodes.in.RWgroup)%>%as.data.frame()%>%mutate(condition="DS"),
  bzip23.tbl%>%select(nodes.in.RWgroup)%>%as.data.frame()%>%mutate(condition="bzip23")
  #bzip23.tbl%>%select(degree.in.RWgroup)%>%as.data.frame()%>%mutate(condition="bzip23")
)

nodes.in.RWgroup.data$nodes.in.RWgroup <- as.numeric(nodes.in.RWgroup.data$nodes.in.RWgroup)
nodes.in.RWgroup.data$condition <- as.factor(nodes.in.RWgroup.data$condition)

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/01_RW_group_size.pdf",bg=bg,height = 4,width = 4)
ggpubr::ggarrange(
ggplot(degree.in.RWgroup.data,aes(x=condition,y=log2(degree.in.RWgroup)))+
  geom_violin(aes(fill=condition),trim = T)+
  geom_boxplot(width=0.2)+
  scale_fill_manual(values = col.conditions)+
  ggsignif::geom_signif(comparisons = list(c("DS","bzip23"))),

ggplot(nodes.in.RWgroup.data,aes(x=condition,y=log2(nodes.in.RWgroup)))+
  geom_violin(aes(fill=condition),trim = T)+
  geom_boxplot(width=0.2)+
  scale_fill_manual(values = col.conditions)+
  ggsignif::geom_signif(comparisons = list(c("DS","bzip23"))),

ncol = 2,align = "hv",common.legend = T
)
dev.off()

#########===============================================================================================#########
#########Characterize promoter-associated super-anchors// method 2
#########===============================================================================================#########

###In this section, we process super-anchors by 3 criteria:
###1. The anchors should also be the hub (with the highest degree, minimal of 3) in it's PPI community (by random walk).
###2. The FG group it belongs to should not be too small, at least the total anchors >= 4 and total degrees >=8 is required.
###3. The hub-anchor should have a least degree of 3.
###4. We define the connectivity density (degree.dens) to evaluate the ability of the anchor to integrate PPIs. hub-anchors with a
###degree.dens over the 75% quantile of all PPI-anchors under the condition are considered as super-anchors.


#Get the hub anchor with the largest degree in the FG group.
#We get all hubs at this step and only filter the degree not group size.
lapply(c("D","bzip23") ,function(x){
  tbl <- get(paste0(x,".tbl"))

  #total.ipet <- tbl%E>%select(ipet.counts)%>%pull()%>%sum()

  tbl%>%as.data.frame()%>%
    group_by(RW.group)%>%
    filter(degree==max(degree,na.rm=T))%>% #get the anchor with the max degree in the RW group.
    ungroup()%>%
    rowwise()%>%
    filter(degree >=3)%>% ##The hub anchor should have a least degree of 3, or it is considered less connected.
    arrange(degree.dens)%>% #Rank the anchors by degree density.
    as.data.frame()%>%
    assign(paste0(x,".HubAnchor"),.,envir = .GlobalEnv)
}
)

D.HubAnchor%>%dim() #722

bzip23.HubAnchor%>%dim() #475


#plot the degree.dens distribution of the hub anchors.
HubAnchor.degree.dens.data <- rbind(
  D.HubAnchor%>%select(anchor.id,degree.dens)%>%mutate(condition="DS"),
  bzip23.HubAnchor%>%select(anchor.id,degree.dens)%>%mutate(condition="bzip23")
)%>%
  group_by(condition)%>%
  arrange(degree.dens)%>%
  mutate(rank=seq(1:n()))%>%
  as.data.frame()%>%
  mutate(rank=as.numeric(rank),condition=as.factor(condition))


#Use the 80% quantile of Hub degree density under the corresponding condition as cutoff

degree.dens.cutoff.D <- quantile(D.HubAnchor$degree.dens,probs = c(0.8))[1]%>%as.vector()
degree.dens.cutoff.bzip23 <- quantile(bzip23.HubAnchor$degree.dens,probs = c(0.8))[1]%>%as.vector()


#########===============================================================================================#########
#########Characterize super-anchors
#########===============================================================================================#########
###Filter the Hub anchors by degree.dens to get super-anchors

lapply(c("D","bzip23") ,function(x){
  hubs <- get(paste0(x,".HubAnchor"))
  #H3K9ac.loaded.anchors <- get(paste0(x,".loaded.anchors"))

  cutoff <- get(paste0("degree.dens.cutoff.",x))

  sa <-hubs%>%
    filter(nodes.in.RWgroup >= 4,degree.in.RWgroup >= 8)%>% ##Discard less connected groups.
    #filter(anchor.id %in% H3K9ac.loaded.anchors)%>% ##Preserve only anchors heavily loaded with H3K9ac.
    filter(degree.dens>cutoff) ##filter contact density

  assign(paste0(x,".superAnchors"),sa,envir = .GlobalEnv)

  write.xlsx(x = sa,file = paste0("../03_Super_Anchor/01_SuperAnchor_list/",x,".superAnchors.xlsx"))

  sa%>%rowwise%>%
    mutate(v1=str_split(anchor.id,pattern = "_",simplify = T)[1],
           v2=str_split(anchor.id,pattern = "_",simplify = T)[2],
           v3=str_split(anchor.id,pattern = "_",simplify = T)[3],
           v4=anchor.id,
           v5=degree.dens,
           v6="*")%>%
    select(v1,v2,v3,v4,v5,v6)%>%
    arrange(v1,v2,v3)%>%
    write.table(paste0("../03_Super_Anchor/02_SuperAnchor_bed/",x,".superAnchors.bed"),
                sep = "\t",row.names = F,col.names = F,quote = F)
}
)

D.superAnchors%>%dim() #130
bzip23.superAnchors%>%dim() #89



#Plot the ranked degree density and trace the positions of the DS super-anchors in NC and RE sets.
#PPIs
lapply(c("D","bzip23"), function(x){

  DSA.mark <- D.superAnchors%>%
    select(anchor.id,degree.dens)%>%
    rename(degree.dens.DS=degree.dens)%>%
    arrange(degree.dens.DS)%>%
    mutate(color.rank=rownames(.))

  df <- get(paste0(x,".tbl"))%>%
    as.data.frame()%>%
    mutate(DS.top=ifelse(anchor.id %in% D.superAnchors$anchor.id,"Yes","NO"))%>%
    arrange(degree.dens)%>%
    mutate(rank=as.numeric(rownames(.)))


  cutoff <- get(paste0("degree.dens.cutoff.",x))

  gg <- ggplot(df,aes(x=rank,y=degree.dens))+
    geom_hline(yintercept = cutoff,color="grey50",linetype=2)+
    geom_point(color="grey",size=1.5)+
    geom_point(data = df%>%filter(DS.top=="Yes")%>%left_join(DSA.mark),
               aes(color=as.numeric(color.rank)),size=1.5)+
    scale_color_continuous(type = "viridis")+
    ylab("Contact density")

  assign(paste0("graph.",x,".degree.dens.rank.allanchors"),gg,envir = .GlobalEnv)

})

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/06_degree_dens_rank.allanchors.pdf",height = 5,width = 4,bg=bg)
ggpubr::ggarrange(graph.D.degree.dens.rank.allanchors,graph.bzip23.degree.dens.rank.allanchors,
                  ncol = 1,align = "hv",common.legend=T)
dev.off()


#Hubs
lapply(c("D","bzip23"), function(x){

  DSA.mark <- D.superAnchors%>%
    select(anchor.id,degree.dens)%>%
    rename(degree.dens.DS=degree.dens)%>%
    arrange(degree.dens.DS)%>%
    mutate(color.rank=rownames(.))

  df <- get(paste0(x,".HubAnchor"))%>%
    as.data.frame()%>%
    mutate(DS.top=ifelse(anchor.id %in% D.superAnchors$anchor.id,"Yes","NO"))%>%
    arrange(degree.dens)%>%
    mutate(rank=as.numeric(rownames(.)))


  cutoff <- get(paste0("degree.dens.cutoff.",x))

  gg <- ggplot(df,aes(x=rank,y=degree.dens))+
    geom_hline(yintercept = cutoff,color="grey50",linetype=2)+
    geom_point(color="grey",size=1.5)+
    geom_point(data = df%>%filter(DS.top=="Yes")%>%left_join(DSA.mark),
               aes(color=as.numeric(color.rank)),size=1.5)+
    scale_color_continuous(type = "viridis")+
    ylab("Contact density")

  assign(paste0("graph.",x,".degree.dens.rank.hub"),gg,envir = .GlobalEnv)

})

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/07_degree_dens_rank.Hubanchors.pdf",height = 5,width = 4,bg=bg)
ggpubr::ggarrange(graph.D.degree.dens.rank.hub,graph.bzip23.degree.dens.rank.hub,
                  ncol = 1,align = "hv",common.legend=T)
dev.off()

#########===============================================================================================#########
#########Identify differentially connected super-anchor.
#########===============================================================================================#########
#1.5-fold change in the degree density value is used to identify a differential connected super-anchor between 2 conditions.

##Use the full quantification list when comparing, the anchor can still get a density value even if it is not a hub anchor under the condition.
##This helps to distinguish true dynamic hubs from constitutively connected hubs.
SuperAnchor.degree.dens.mtx <- rbind(D.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="DS"),
                                     bzip23.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="bzip23"))%>%
  filter(anchor.id %in% c(D.superAnchors$anchor.id,bzip23.superAnchors$anchor.id))%>%
  tidyr::pivot_wider(id_cols=c(anchor.id),
                     names_from =c(condition),
                     values_from =c(degree.dens),
                     values_fill =list(degree.dens=0.001))%>% #give a pseudo density to avoid error if the density is absent under the condition.
  as.data.frame()


#Get the DS super anchors who became less connected in bzip23.
#The dominate anchor here is in terms of DS vs bzip23.
bzip23.affected.DS.SuperAnchor <- SuperAnchor.degree.dens.mtx %>%
  filter(DS/bzip23>=1.5)%>%
  select(-DS,-bzip23)%>%
  filter(anchor.id %in% D.superAnchors$anchor.id)%>%
  left_join(PPI.anchor.info) #58 vs 130

#See if the mutation of bZIP23 can further improve the chromation loop-integrating ability of DS super-anchors
SuperAnchor.degree.dens.mtx %>%
  filter(bzip23/DS>=1.5)%>%
  select(-DS,-bzip23)%>%
  filter(anchor.id %in% D.superAnchors$anchor.id)%>%dim()# 3 vs 130
#The answer is no.

#See how many bzip23 super-anchors are specific to bzip23 rather than DS.
bzip23.denovo.SuperAnchor <- SuperAnchor.degree.dens.mtx%>%
  filter(bzip23/DS>=1.5)%>%
  select(-DS,-bzip23)%>%
  filter(anchor.id %in% bzip23.superAnchors$anchor.id)%>%
  left_join(PPI.anchor.info) #53 vs 89


#See how many DS-dominate super anchors are affected by bzip23 mutation.
D.condition.domi.superAnchor <- read.xlsx("../../01_analysis_on_MH63WT/03_Super_Anchor/03_Dominate_superAnchor_list/DS_dominate_superAnchor.xlsx")

bzip23.affected.DS.condi.domi.superAnchor <- SuperAnchor.degree.dens.mtx%>%
  filter(DS/bzip23>=1.5)%>%
  select(-DS,-bzip23)%>%
  filter(anchor.id %in% D.condition.domi.superAnchor$anchor.id)%>%
  left_join(PPI.anchor.info) #37 vs 68

#Next see how many bZIP23-affected DS-dominate super-anchors are bound directly by bZIP23.
bzip23.bind.anchor <- read.xlsx("../../../25_bzip23_TFBS_MH63_ChIApet_integrative/bZIP23_bind_anchor_DS.xlsx")%>%
  select(V4,V5)%>%
  rename(anchor.id=V4,bZIP23.bind.count=V5)


bzip23.affected.DS.SuperAnchor%<>%left_join(bzip23.bind.anchor)
bzip23.affected.DS.SuperAnchor$bZIP23.bind.count[is.na(bzip23.affected.DS.SuperAnchor$bZIP23.bind.count)]<-0

bzip23.affected.DS.SuperAnchor%>%filter(bZIP23.bind.count>0)%>%dim() #28 vs 58
dhyper(28,5243,21065-5243,58) #6.511492e-05

bzip23.affected.DS.condi.domi.superAnchor%<>%left_join(bzip23.bind.anchor)
bzip23.affected.DS.condi.domi.superAnchor$bZIP23.bind.count[is.na(bzip23.affected.DS.condi.domi.superAnchor$bZIP23.bind.count)]<-0

bzip23.affected.DS.condi.domi.superAnchor%>%filter(bZIP23.bind.count>0)%>%dim() #18 vs 37 #0.4864865#==================
dhyper(18,5243,21065-5243,37)  #0.001023968

#The bzip23 affected super-anchors bound by bZIP23 are regarded as bZIP23-regulating super-anchors.
bzip23.reg.DS.SuperAnchor <- bzip23.affected.DS.SuperAnchor%>%filter(bZIP23.bind.count>0)

bzip23.reg.DS.condi.domi.SuperAnchor <- bzip23.affected.DS.condi.domi.superAnchor%>%filter(bZIP23.bind.count>0)

##Also see how many DS-PPI-anchors and DS-super-anchors are bound by bZIP23
#DS-PPI-anchors
intersect(D.tbl%>%as.data.frame()%>%pull(anchor.id),bzip23.bind.anchor$anchor.id)%>%length()# 2396 / 6091 #0.3933673

#DS-super-anchors
D.superAnchors%>%left_join(bzip23.bind.anchor)%>%filter(is.na(bZIP23.bind.count)==F)%>%dim()# 70 / 130 #0.5384615

#DS-dominate-super-anchor
SuperAnchor.degree.dens.mtx%>%
  #filter(DS/bzip23>=1.5)%>%
  filter(anchor.id %in% D.condition.domi.superAnchor$anchor.id)%>% #68
  left_join(bzip23.bind.anchor)%>%filter(is.na(bZIP23.bind.count)==F)%>%dim() #33 / 68 #0.4852941#=================

####May 17 2024 (GB rev.1):Find the correlation of DS/bZIP23 CIA value and bZIP23 binding strength (peak intensities).
###Load bZIP23 ChIP-seq quantification matrix on DS-specific SPRs
load("bZIP23_quant_on_bed.RData")

bZIP23_ChIP_quant_on_DS_dominate_SPR <- mtx.count$binding%>%as.data.frame()%>%
                                        cbind(.,D.condition.domi.superAnchor%>%select(anchor.id))%>%
                                        mutate(NC=(NC_1+NC_2)/2,DS=(DS_1+DS_2)/2,ratio=DS/NC)

df1 <- bZIP23_ChIP_quant_on_DS_dominate_SPR%>%filter(anchor.id %in% bzip23.affected.DS.condi.domi.superAnchor$anchor.id)
df2 <- bZIP23_ChIP_quant_on_DS_dominate_SPR%>%filter(anchor.id %in% D.condition.domi.superAnchor$anchor.id)

wilcox.test(df1$DS,df2$DS)

pdf("./bZIP23_peak_intens_on_DS_dominate_SPR.pdf",height = 4,width = 3,bg=bg)
rbind(data.frame(value=df1$DS,cat="DS_specific_SPR_lost_in_osbzip23"),data.frame(value=df2$DS,cat="DS_specific_SPR"))%>%
ggplot(aes(x=cat,y=value))+
  geom_boxplot(outlier.shape = NA)+
  geom_signif(comparisons = list(c("DS_specific_SPR","DS_specific_SPR_lost_in_osbzip23")))+
  xlab("")+
  ylab("log2 OsbZIP23 binding intensity (IP vs input)")
dev.off()

##plot heatmap of the degree.dens matrix
library(viridis)

heatmap.in <- SuperAnchor.degree.dens.mtx%>%
  mutate(is.DS.condi.domi.SA = ifelse(anchor.id %in% D.condition.domi.superAnchor$anchor.id,"Yes","NO"))%>%
  mutate(is.bzip23.affected.DS.condi.domi.SA = ifelse(anchor.id %in% bzip23.affected.DS.condi.domi.superAnchor$anchor.id,"Yes","NO"))%>%
  mutate(is.bzip23.reg.DS.condi.domi.SA = ifelse(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor$anchor.id,"Yes","NO"))

rownames(heatmap.in) <- heatmap.in$anchor.id

pdf(file = "DS_vs_bzip23.AnchorDense.heatmap.pdf",height = 7,width =7 ,bg = bg)
pheatmap::pheatmap(
  heatmap.in%>%select(-anchor.id,-is.DS.condi.domi.SA,-is.bzip23.affected.DS.condi.domi.SA,-is.bzip23.reg.DS.condi.domi.SA),
  col=viridis::inferno(n = 255,alpha = 1,direction = -1),
  annotation_row=heatmap.in%>%select(is.DS.condi.domi.SA,is.bzip23.affected.DS.condi.domi.SA,is.bzip23.reg.DS.condi.domi.SA),
  show_rownames = F,
  clustering_method = "complete",
  annotation_colors  = list(
    is.DS.condi.domi.SA=c(Yes=col.DS,NO="white"),
    is.bzip23.affected.DS.condi.domi.SA=c(Yes=col.bzip23,NO="white"),
    is.bzip23.reg.DS.condi.domi.SA=c(Yes=col.bzip23,NO="white"))
)
dev.off()

##Export theses anchors
#Spread sheet
write.xlsx(bzip23.affected.DS.SuperAnchor,"../03_Super_Anchor/01_SuperAnchor_list/bzip23.affected.DS.SuperAnchor.xlsx")
write.xlsx(bzip23.affected.DS.condi.domi.superAnchor,"../03_Super_Anchor/01_SuperAnchor_list/bzip23.affected.DS_condition_dominate.SuperAnchor.xlsx")
write.xlsx(bzip23.reg.DS.condi.domi.SuperAnchor,"../03_Super_Anchor/01_SuperAnchor_list/bzip23.regulated.DS_condition_dominate.SuperAnchor.xlsx")


#########===============================================================================================#########
#########Analyze the effect of bzip23-regulated SPE in the integration of RW group
#########===============================================================================================#########
####Compare the log2FC of PPI genes by DS-dominate PPIs integrated by bzip23-reg dominate SPE and non-bzip23-reg ones.
group.wi.bzip23.reg.SPE <- D.tbl%>%filter(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor$anchor.id)%>%
  as.data.frame()%>%
  select(RW.group)

group.wi.ds.domi.PSE <- D.tbl%>%filter(anchor.id %in% D.condition.domi.superAnchor$anchor.id)%>%
  as.data.frame()%>%
  select(RW.group)


pdf("../01_PPI_network_analysis/02_PPI_network_profiles/08_anchor_logFC_by_bzip23_affSPE_group.pdf",height = 4,width = 4,bg=bg)
rbind(
  D.tbl%>%
    filter(RW.group %in% group.wi.bzip23.reg.SPE$RW.group)%E>%
    filter(loop.domination=="DS_dominate")%N>%
    filter(!(node_is_isolated()))%N>%
    select(DS.TPM.anchor,avgTPM.DS.RWgroup,avgTPM.NC.RWgroup,avgTPM.RE.RWgroup,log2FC_DS_vs_NC.anchor)%>%
    as.data.frame()%>%
    mutate(cat="bzip23.reg.SPE.RWgroup"),

  D.tbl%>%
    filter(!(RW.group %in% group.wi.bzip23.reg.SPE$RW.group),
           RW.group %in%group.wi.ds.domi.PSE$RW.group)%E>%
    filter(loop.domination=="DS_dominate")%N>%
    filter(!(node_is_isolated()))%N>%
    select(DS.TPM.anchor,avgTPM.DS.RWgroup,avgTPM.NC.RWgroup,avgTPM.RE.RWgroup,log2FC_DS_vs_NC.anchor)%>%
    as.data.frame()%>%
    mutate(cat="non-bzip23.reg.SPE.RWgroup")
)%>%
  ggplot(aes(x=cat,y=log2FC_DS_vs_NC.anchor))+
  geom_violin(draw_quantiles = F,trim = F)+
  geom_boxplot(width=0.3,outlier.shape = NA)+
  geom_signif(comparisons = list(c("bzip23.reg.SPE.RWgroup","non-bzip23.reg.SPE.RWgroup")))
dev.off()

######Compare the contact intensities between bzip23 reg dominate SPEs and non-bzip23 reg ones.

  rbind(

    D.HubAnchor%>%filter(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor$anchor.id)%>%
      select(degree.dens,anchor.id)%>%
      mutate(anchor.type="bzip23 reg dominate SPE")%>%unique(),

    D.HubAnchor%>%filter(anchor.id %in% D.condition.domi.superAnchor$anchor.id)%>%
      filter(!(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor$anchor.id))%>%
      select(degree.dens,anchor.id)%>%
      mutate(anchor.type="other dominate SPE")%>%unique()
  )%>%
  ggplot(aes(x=anchor.type,y=degree.dens))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  geom_signif(comparisons = list(c("bzip23 reg dominate SPE","other dominate SPE")))

  D.tbl%N>%filter(RW.group %in% group.wi.bzip23.reg.SPE$RW.group)%>%
    as.data.frame()%>%select(name,RW.group,degree.dens)%>%
    unique()%>%filter(name != "")%>%
    left_join(D.tbl%N>%filter(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor$anchor.id)%>%as.data.frame()%>%select(anchor.id,RW.group))%>%na.omit

#########===============================================================================================#########
#########Plot PPI network of DS and bzip23.
#########===============================================================================================#########
node.column <- c("anchor.id","degree","degree.dens","represent.gene","DS.TPM.anchor","NC.TPM.anchor",
                 "RE.TPM.anchor","bzip23.TPM.anchor","log2FC_DS_vs_NC.anchor","log2FC_RE_vs_DS.anchor",
                 "log2FC_RE_vs_NC.anchor","log2FC_DS_vs_bzip23.anchor",
                 "mfuzz.cluster.anchor","membership.mfuzz.1","membership.mfuzz.2","membership.mfuzz.3",
                 "membership.mfuzz.4","membership.mfuzz.5","membership.mfuzz.6","MSU7_ID","name")

anchor.label <- c("OsjDHN6|RAB16C","OsjDHN7|RAB16B","OsjDHN5","OsjDHN6|RAB16C")

tbl <- graph_join(D.tbl%N>%select(all_of(node.column))%>%rename(degree.DS=degree,degree.dens.DS=degree.dens),
                  bzip23.tbl%N>%select(all_of(node.column))%>%rename(degree.bzip23=degree,degree.dens.bzip23=degree.dens))

tbl %<>% mutate(is.DS.domi.superAnchor = ifelse(anchor.id %in% D.condition.domi.superAnchor$anchor.id,"Yes","NO"),
                is.bzip23.reg.DS.domi.SPE = ifelse(anchor.id %in% bzip23.reg.DS.condi.domi.SuperAnchor,"Yes","NO"))

tbl.chr11 <- tbl%>%filter(str_split(anchor.id,"_",simplify = T)[,1]=="Chr11")%>%filter(!(node_is_isolated()))

set.seed(100)
tbl.chr11.weights <-tbl.chr11%E>%as.data.frame%>%select(ipet.counts)%>%pull()
tbl.chr11.xy <- layout_with_fr(tbl.chr11,grid="nogrid")


#--plot DS
pdf("../01_PPI_network_analysis/03_PPI_network_plot/DS.allPPIs.Chr11.pdf",height = 15,width = 15,bg=bg)
ggraph(tbl.chr11,layout = "manual",x=tbl.chr11.xy[,1],y=tbl.chr11.xy[,2])+
  geom_edge_density(aes(fill=loop.domination,
                        filter=loop.in.DS==1))+
  scale_edge_fill_manual(values=c("other"=bg,
                                  "Conserved"=bg,
                                  "RE_dominate"=col.RE,
                                  "NC_dominate"=col.NC,
                                  "DS_dominate"=col.DS))+
  geom_edge_link(aes(edge_color=loop.domination,
                     edge_width=loop.domination,
                     filter=loop.in.DS==1),
                 check_overlap =T)+
  scale_edge_width_manual(values=c("other"=0.8,
                                   "Conserved"=0.8,
                                   "RE_dominate"=2,
                                   "NC_dominate"=2,
                                   "DS_dominate"=2))+
  scale_edge_color_manual(values=c("other"="grey75",
                                   "Conserved"="grey70",
                                   "RE_dominate"=bg,
                                   "NC_dominate"=bg,
                                   "DS_dominate"=col.DS))+
  geom_node_point(aes(fill=as.factor(mfuzz.cluster.anchor),
                      size=membership.mfuzz.2),
                  shape=21)+
  scale_size_continuous(range = c(2,8))+
  #highlight DS-dominate super-anchors
  geom_node_point(aes(filter=is.DS.domi.superAnchor=="Yes"),
                  shape=23,size=8,fill=col.DS)+
  geom_node_point(aes(filter=is.bzip23.reg.DS.domi.SPE=="Yes"),
                  shape=22,size=8,fill=col.DS)+
  scale_fill_manual(values=c("1"="lightgreen",
                             "2"="pink",
                             "3"="pink",
                             "4"="cyan",
                             "5"="lightgreen",
                             "6"="skyblue",
                             "NA"="grey50"))+
  geom_node_text(aes(label=name,filter=name %in% anchor.label),size=5)+
  theme_graph(base_family= "Arial")+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(title="")
dev.off()


#--plot bzip23
pdf("../01_PPI_network_analysis/03_PPI_network_plot/bzip23.allPPIs.Chr11.pdf",height = 15,width = 15,bg=bg)
ggraph(tbl.chr11,layout = "manual",x=tbl.chr11.xy[,1],y=tbl.chr11.xy[,2])+
  geom_edge_density(aes(fill=loop.domination,
                        filter=loop.in.bzip23==1))+
  scale_edge_fill_manual(values=c("other"=bg,
                                  "Conserved"=bg,
                                  "RE_dominate"=col.RE,
                                  "NC_dominate"=col.NC,
                                  "DS_dominate"=col.DS))+
  geom_edge_link(aes(edge_color=loop.domination,
                     edge_width=loop.domination,
                     filter=loop.in.bzip23==1),
                 check_overlap =T)+
  scale_edge_width_manual(values=c("other"=0.8,
                                   "Conserved"=0.8,
                                   "RE_dominate"=2,
                                   "NC_dominate"=2,
                                   "DS_dominate"=2))+
  scale_edge_color_manual(values=c("other"="grey75",
                                   "Conserved"="grey70",
                                   "RE_dominate"=col.RE,
                                   "NC_dominate"=col.NC,
                                   "DS_dominate"=col.DS))+
  geom_node_point(aes(fill=as.factor(mfuzz.cluster.anchor),
                      size=membership.mfuzz.2),
                  shape=21)+
  scale_size_continuous(range = c(2,8))+
  #highlight DS-dominate super-anchors
  geom_node_point(aes(filter=is.DS.domi.superAnchor=="Yes"),
                  shape=23,size=8,fill=col.DS)+
  geom_node_point(aes(filter=is.bzip23.reg.DS.domi.SPE=="Yes"),
                  shape=22,size=8,fill=col.DS)+
  scale_fill_manual(values=c("1"="lightgreen",
                             "2"="pink",
                             "3"="pink",
                             "4"="cyan",
                             "5"="lightgreen",
                             "6"="skyblue",
                             "NA"="grey50"))+
  geom_node_text(aes(label=name,filter=name %in% anchor.label),size=5)+
  theme_graph(base_family= "Arial")+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(title="")
dev.off()
