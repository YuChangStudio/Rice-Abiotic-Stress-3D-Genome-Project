#########===============================================================================================#########
#########In this section, we analyze putative enhancer/super-promoter associating anchors by network analysis
#########===============================================================================================#########

#########===============================================================================================#########
#########Load configs and data
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

col.conditions <- c("DS"=col.DS,"NC"=col.NC,"RE"=col.RE)

#ggplot2 theme

theme(
      panel.background = element_rect(fill = "transparent",colour='black',size = 1),
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


##Load only the PPI loops, other loops are not applicable for this analysis.
##Use the PPI data with unique representative gene info, since a node (anchor) can be labeled with only one set of info.
##Inter-chromosomal interactions are not considered in the analysis,
##This unique gene annotation files were generated
##The PPI object will be used as edge input.
##The anchor annotation file with representative gene will be used as node input.

##Generate network node object:
PPI.anchor.info <- read.xlsx("../../../20_ChIAPET_profile_analysis/03_PPI_pairs/03_PPI_anchor_annoation/PPI_anchor_unique_genes.annoation.xlsx")%>%
                left_join(.,mh63.encode,by=c("represent.gene"="MH63RS3_ID"))%>%
                mutate(node.name=ifelse(name != "" , name, represent.gene))%>%
                select(-bzip23.TPM.anchor,-log2FC_DS_vs_bzip23.anchor)

##Load results on condition-dominate PPIs as edge info.
PPI.domination.category <- read.xlsx("../../../22_diffloop_analysis/01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")%>%
                            mutate(loop.domination=case_when(
                              (loop.in.DS ==1 & loop.in.NC==0 & loop.in.RE==0) ~ "DS_dominate",
                              (loop.in.DS ==0 & loop.in.NC==1 & loop.in.RE==0) ~ "NC_dominate",
                              (loop.in.DS ==0 & loop.in.NC==0 & loop.in.RE==1) ~ "RE_dominate",
                              (loop.in.DS ==1 & loop.in.NC==1 & loop.in.RE==1) ~ "Conserved",
                              TRUE ~ "other"
                            ))%>%
                            select(-anchor1.id,-anchor2.id)

#########===============================================================================================#########
#########Generate tidygraph network object and make basic analysis
#########===============================================================================================#########
##Generate tbl network object:
#degree.dens=(degree * 1000 * 10000)/(total degree in data set * anchor length in bp)
#We expect that for more connected data set, the super anchor should have larger degree.
lapply(c("D","N","RE"), function(x){

  file.prefix <- "../../../20_ChIAPET_profile_analysis/03_PPI_pairs/04_PPI_loops_with_representative_genes/"

  #get raw edge
  PPI <- read.xlsx(paste0(file.prefix,x,".PPI.uniqGene.xlsx"))%>%
    select(anchor1.id,anchor2.id,loop.index,ipet.counts)%>%#get edge info
    left_join(.,PPI.domination.category,by="loop.index")

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
           medianTPM.RE.RWgroup=median(RE.TPM.anchor)
           )%>%
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

##FG/LP/RW group number under each condition
#NC.tbl%>%select(FG.group)%>%as.data.frame()%>%max()#1005
#NC.tbl%>%select(PL.group)%>%as.data.frame()%>%max()#2060
N.tbl%>%select(RW.group)%>%as.data.frame()%>%max()#2097

#DS.tbl%>%select(FG.group)%>%as.data.frame()%>%max()#1046
#DS.tbl%>%select(PL.group)%>%as.data.frame()%>%max()#1630
D.tbl%>%select(RW.group)%>%as.data.frame()%>%max()#1427

#RE.tbl%>%select(FG.group)%>%as.data.frame()%>%max()#976
#RE.tbl%>%select(PL.group)%>%as.data.frame()%>%max()#1755
RE.tbl%>%select(RW.group)%>%as.data.frame()%>%max()#1674

####Compare the RW group size (by degree and node numbers):
RWgroup.size.data <- rbind(
N.tbl%>%select(degree.in.RWgroup,nodes.in.RWgroup)%>%as.data.frame()%>%mutate(condition="NC"),
D.tbl%>%select(degree.in.RWgroup,nodes.in.RWgroup)%>%as.data.frame()%>%mutate(condition="DS"),
RE.tbl%>%select(degree.in.RWgroup,nodes.in.RWgroup)%>%as.data.frame()%>%mutate(condition="RE")
)

#degree.in.RWgroup.data$degree.in.RWgroup <- as.numeric(degree.in.RWgroup.data$degree.in.RWgroup)
#degree.in.RWgroup.data$condition <- as.factor(degree.in.RWgroup.data$condition)
#HSD test
#degree.in.RWgroup.hsd <- lm(log2(degree.in.RWgroup)~condition,data = degree.in.RWgroup.data)%>%
#  multcomp::glht(.,linfct = multcomp::mcp(condition = "Tukey"))%>%multcomp::cld(decreasing=T)

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/01_RWgroup_size.violin.pdf",bg=bg,height = 4,width = 6)
ggplot(RWgroup.size.data%>%reshape2::melt(),aes(x=condition,y=log2(value)))+
  geom_violin(aes(fill=condition),trim = F)+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values = col.conditions)+
  geom_signif(comparisons = list(c("NC","DS"),c("NC","RE"),c("RE","DS")),step_increase = 0.1)+
  facet_wrap(.~variable,ncol = 2,scales = "free")
dev.off()

##node vs degree
##exclude the largest group since it will affect the linear fit greatly.
pdf("../01_PPI_network_analysis/02_PPI_network_profiles/01_RWgroup_node_vs_degree.scatter.pdf",bg=bg,height = 4,width = 6)
ggplot(RWgroup.size.data%>%filter(degree.in.RWgroup<1000),aes(x=nodes.in.RWgroup,y=degree.in.RWgroup))+
  geom_point(aes(color=condition))+
  geom_smooth(aes(group=condition),method = "lm")+
  scale_color_manual(values = col.conditions)+
  facet_wrap(.~condition,ncol = 3,scales = "fixed")
dev.off()

lm(degree.in.RWgroup~nodes.in.RWgroup,data = RWgroup.size.data%>%filter(degree.in.RWgroup<1000)%>%filter(condition=="NC")) #y=9.05*x-43.40
lm(degree.in.RWgroup~nodes.in.RWgroup,data = RWgroup.size.data%>%filter(degree.in.RWgroup<1000)%>%filter(condition=="DS")) #y=4.707*x-11.693
lm(degree.in.RWgroup~nodes.in.RWgroup,data = RWgroup.size.data%>%filter(degree.in.RWgroup<1000)%>%filter(condition=="RE")) #y=5.612*x-11.554

##Plot the correlation of degree and anchor span.

degree_vs_range.data <- rbind(
  N.tbl%>%select(degree,anchor.range)%>%as.data.frame()%>%mutate(condition="NC"),
  D.tbl%>%select(degree,anchor.range)%>%as.data.frame()%>%filter(degree < 120)%>%mutate(condition="DS"),
  RE.tbl%>%select(degree,anchor.range)%>%as.data.frame()%>%mutate(condition="RE")
)

cor(N.tbl%>%select(degree,anchor.range)%>%as.data.frame()) # 0.8053981
cor(D.tbl%>%select(degree,anchor.range)%>%as.data.frame()%>%filter(degree < 120)) # 0.735366
cor(RE.tbl%>%select(degree,anchor.range)%>%as.data.frame()) # 0.7694882

degree_vs_range.graph<-ggplot(degree_vs_range.data,aes(x=degree,y=anchor.range))+
  geom_point(aes(color=condition))+
  geom_smooth(aes(group=condition),method = "lm")+
  scale_color_manual(values = col.conditions)+
  facet_wrap(.~condition,ncol = 3,scales = "free")

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/02_degree_vs_anchor_span.pdf",bg=bg,height = 4,width = 10)
degree_vs_range.graph
dev.off()

##plot the correlation of degree.dens vs anchor span
degree.dens_vs_range.data <- rbind(
  N.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()%>%mutate(condition="NC"),
  D.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()%>%mutate(condition="DS"),
  RE.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()%>%mutate(condition="RE")
)

cor(N.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()) # 0.185878
cor(D.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()) # -0.2550582
cor(RE.tbl%>%select(degree.dens,anchor.range)%>%as.data.frame()) # -0.02758044

degree.dens_vs_range.graph <- ggplot(degree.dens_vs_range.data,aes(x=degree.dens,y=anchor.range))+
  geom_point(aes(color=condition))+
  geom_smooth(aes(group=condition),method = "lm")+
  scale_color_manual(values = col.conditions)+
  facet_wrap(.~condition,ncol = 3,scales = "free")

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/03_degree.dens_vs_anchor_span.pdf",bg=bg,height = 4,width = 10)
degree.dens_vs_range.graph
dev.off()

##plot the correlation of degree.dens vs degree
degree.dens_vs_degree.data <- rbind(
  N.tbl%>%select(degree.dens,degree)%>%as.data.frame()%>%mutate(condition="NC"),
  D.tbl%>%select(degree.dens,degree)%>%as.data.frame()%>%mutate(condition="DS"),
  RE.tbl%>%select(degree.dens,degree)%>%as.data.frame()%>%mutate(condition="RE")
)

cor(N.tbl%>%select(degree.dens,degree)%>%as.data.frame()) # 0.555506
cor(D.tbl%>%select(degree.dens,degree)%>%as.data.frame()) # 0.2919742
cor(RE.tbl%>%select(degree.dens,degree)%>%as.data.frame()) # 0.436797


degree.dens_vs_degree.graph <- ggplot(degree.dens_vs_degree.data,aes(x=degree.dens,y=degree))+
  geom_point(aes(color=condition))+
  geom_smooth(aes(group=condition),method = "lm")+
  scale_color_manual(values = col.conditions)+
  facet_wrap(.~condition,ncol = 3,scales = "free")

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/04_degree.dens_vs_degree.pdf",bg=bg,height = 4,width = 10)
degree.dens_vs_degree.graph
dev.off()

##Combine the three plots
pdf("../01_PPI_network_analysis/02_PPI_network_profiles/05_degree.dens_vs_degree_vs_span.pdf",bg=bg,height = 10,width = 10)
ggpubr::ggarrange(degree_vs_range.graph,
                  degree.dens_vs_range.graph,
                  degree.dens_vs_degree.graph,
                  ncol = 1,
                  nrow = 3,
                  align = "hv",
                  legend = "right")
dev.off()

#########===============================================================================================#########
#########Characterize promoter-associated super-anchors// method 2
#########===============================================================================================#########

###In this section, we process super-anchors by 3 criteria:
###1. The anchors should be the hub (with the highest degree, minimal of 3) in it's PPI community (by random walk).
###2. The RW group it belongs to should not be too small, at least the total anchors >= 4 and total degrees >=8 are required.
###3. The hub-anchor should have a least degree of 3.
###4. We define the connectivity density (degree.dens) to evaluate the ability of the anchor to integrate PPIs. hub-anchors with a
###degree.dens over the 75% quantile of all PPI-anchors under the condition are considered as super-anchors.


#########===============================================================================================#########
#########Characterize hub PPI anchors
#########===============================================================================================#########

#Get the hub anchor with the largest degree in the RW group.
#We get all hubs at this step and only filter the degree not group size.
lapply(c("D","N","RE") ,function(x){
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
N.HubAnchor%>%dim() #1362
RE.HubAnchor%>%dim() #830

#Chr11_17997217_18000215

#Get the degree.dens distribution of all PPI anchors and the hub anchors, respectively.
PPIAnchor.degree.dens.data <- rbind(
  D.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="DS"),
  N.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="NC"),
  RE.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="RE")
)%>%
  group_by(condition)%>%
  arrange(degree.dens)%>%
  mutate(rank=seq(1:n()))%>%
  as.data.frame()%>%
  mutate(rank=as.numeric(rank),condition=as.factor(condition))

HubAnchor.degree.dens.data <- rbind(
    D.HubAnchor%>%select(anchor.id,degree.dens)%>%mutate(condition="DS"),
    N.HubAnchor%>%select(anchor.id,degree.dens)%>%mutate(condition="NC"),
    RE.HubAnchor%>%select(anchor.id,degree.dens)%>%mutate(condition="RE")
)%>%
group_by(condition)%>%
  arrange(degree.dens)%>%
  mutate(rank=seq(1:n()))%>%
  as.data.frame()%>%
  mutate(rank=as.numeric(rank),condition=as.factor(condition))

quantile(N.HubAnchor%>%select(degree.dens)%>%pull(),probs = c(0.8,0.85,0.9))
quantile(D.HubAnchor%>%select(degree.dens)%>%pull(),probs = c(0.8,0.85,0.9))
quantile(RE.HubAnchor%>%select(degree.dens)%>%pull(),probs = c(0.8,0.85,0.9))


##Use the 80% quantile of Hub degree density under the corresponding condition as cutoff
degree.dens.cutoff.N <- quantile(N.HubAnchor$degree.dens,probs = c(0.8))[1]%>%as.vector()
degree.dens.cutoff.D <- quantile(D.HubAnchor$degree.dens,probs = c(0.8))[1]%>%as.vector()
degree.dens.cutoff.RE <- quantile(RE.HubAnchor$degree.dens,probs = c(0.8))[1]%>%as.vector()

#########===============================================================================================#########
#########Characterize super-anchors
#########===============================================================================================#########
###Filter the Hub anchors by degree.dens to get super-anchors

lapply(c("D","N","RE") ,function(x){
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
N.superAnchors%>%dim() #231
RE.superAnchors%>%dim() #121


###Plot the ranked degree density and trace the positions of the DS super-anchors in NC and RE sets.
##involving all PPIs as background
lapply(c("D","N","RE"), function(x){

  df <- get(paste0(x,".tbl"))%>%
    as.data.frame()%>%
    mutate(DS.top=ifelse(anchor.id %in% D.superAnchors$anchor.id,"Yes","NO"))%>%
    arrange(degree.dens)%>%
    mutate(rank=as.numeric(rownames(.)))


  cutoff <- get(paste0("degree.dens.cutoff.",x))

  gg <- ggplot(df,aes(x=rank,y=degree.dens))+
        geom_hline(yintercept = cutoff,color="grey50",linetype=2)+
        geom_point(color="grey",size=0.1)+
        geom_point(data = df%>%filter(DS.top=="Yes"),color=col.DS,size=0.1)+
        ylab("Contact density")

  assign(paste0("graph.",x,".degree.dens.rank"),gg,envir = .GlobalEnv)

})

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/06_degree_dens_rank.PPIanchors.pdf",height = 6,width = 4,bg=bg)
ggpubr::ggarrange(graph.N.degree.dens.rank,graph.D.degree.dens.rank,graph.RE.degree.dens.rank,
                  ncol = 1,align = "hv")
dev.off()

##involving only Hubs as background
lapply(c("D","N","RE"), function(x){

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

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/07_degree_dens_rank.Hubanchors.pdf",height = 6.2,width = 4,bg=bg)
ggpubr::ggarrange(graph.N.degree.dens.rank.hub,graph.D.degree.dens.rank.hub,graph.RE.degree.dens.rank.hub,
                  ncol = 1,align = "hv",common.legend=T)
dev.off()

#########===============================================================================================#########
#########Identify differentially connected super-anchor.
#########===============================================================================================#########
#1.5-fold change in the degree density value is used to identify a differential connected super-anchor between 2 conditions.


##Use the full quantification list when comparing, the anchor can still get a density value even if it is not a hub anchor under the condition.
##This helps to distinguish true dynamic hubs from constitutively connected hubs.
SuperAnchor.degree.dens.mtx <- rbind(D.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="DS"),
                                     N.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="NC"),
                                     RE.tbl%>%as.data.frame()%>%select(anchor.id,degree.dens)%>%mutate(condition="RE"))%>%
  filter(anchor.id %in% c(D.superAnchors$anchor.id,N.superAnchors$anchor.id,RE.superAnchors$anchor.id))%>%
  tidyr::pivot_wider(id_cols=c(anchor.id),
                     names_from =c(condition),
                     values_from =c(degree.dens),
                     values_fill =list(degree.dens=0.001))%>% #give a pseudo density to avoid error if the density is absent under the condition.
  as.data.frame()


#Get the DS dominate super-anchors (bzip23 not involved in comparison)
D.domi.SuperAnchor <- SuperAnchor.degree.dens.mtx %>%
  filter(DS/NC>=1.5,DS/RE>=1.5)%>%
  select(-DS,-NC,-RE)%>%
  filter(anchor.id %in% D.superAnchors$anchor.id)%>%
  left_join(PPI.anchor.info) #68 vs 130

#Get the NC dominate super-anchors (bzip23 not involved in comparison)
N.domi.SuperAnchor <- SuperAnchor.degree.dens.mtx%>%
  filter(NC/DS>=1.5,NC/RE>=1.5)%>%
  select(-DS,-NC,-RE)%>%
  filter(anchor.id %in% N.superAnchors$anchor.id)%>%
  left_join(PPI.anchor.info) #55 vs 231

#Get the RE dominate super-anchors (bzip23 not involved in comparison)
RE.domi.SuperAnchor <- SuperAnchor.degree.dens.mtx%>%
  filter(RE/DS>=1.5,RE/NC>=1.5)%>%
  select(-DS,-NC,-RE)%>%
  filter(anchor.id %in% RE.superAnchors$anchor.id)%>%
  left_join(PPI.anchor.info) #35 vs 121

##plot heatmap of the degree.dens matrix
library(viridis)
heatmap.in <- SuperAnchor.degree.dens.mtx%>%
  mutate(is.DS.Dominate.SA = ifelse(anchor.id %in% D.domi.SuperAnchor$anchor.id,"Yes","NO"))%>%
  mutate(is.NC.Dominate.SA = ifelse(anchor.id %in% N.domi.SuperAnchor$anchor.id,"Yes","NO"))%>%
  mutate(is.RE.Dominate.SA = ifelse(anchor.id %in% RE.domi.SuperAnchor$anchor.id,"Yes","NO"))

rownames(heatmap.in) <- heatmap.in$anchor.id

pdf(file = "../03_Super_Anchor/05_SuperAnchor_compare_plot/AnchorDense.heatmap.wi.rownames.pdf",height = 7,width =7 ,bg = bg)
pheatmap::pheatmap(
  heatmap.in%>%select(-anchor.id,-is.DS.Dominate.SA,-is.NC.Dominate.SA,-is.RE.Dominate.SA),
  col=viridis::inferno(n = 255,alpha = 1,direction = -1),
  annotation_row=heatmap.in%>%select(is.DS.Dominate.SA,is.NC.Dominate.SA,is.RE.Dominate.SA),
  show_rownames = T,
  clustering_method = "complete",
  annotation_colors  = list(
    is.DS.Dominate.SA=c(Yes=col.DS,NO="white"),
    is.NC.Dominate.SA=c(Yes=col.NC,NO="white"),
    is.RE.Dominate.SA=c(Yes=col.RE,NO="white"))
  )
dev.off()

#Get distance of clustering on samples:
sample.hc <- heatmap.in%>%select(DS,NC,RE)%>%t()%>%dist()%>%hclust(method = "complete")

#Get correlations of the samples based on the degree density.
pdf(file = "../03_Super_Anchor/05_SuperAnchor_compare_plot/AnchorDense.PCCheatmap.pdf",height = 3,width =3 ,bg = bg)
heatmap.in%>%select(DS,NC,RE)%>%cor()%>%pheatmap::pheatmap(display_numbers = T)
dev.off()

##Export the condition-dominate super-anchors
write.xlsx(D.domi.SuperAnchor,"../03_Super_Anchor/03_Dominate_superAnchor_list/DS_dominate_superAnchor.xlsx")
write.xlsx(N.domi.SuperAnchor,"../03_Super_Anchor/03_Dominate_superAnchor_list/NC_dominate_superAnchor.xlsx")
write.xlsx(RE.domi.SuperAnchor,"../03_Super_Anchor/03_Dominate_superAnchor_list/RE_dominate_superAnchor.xlsx")

##Export BedGraphs with 500 bp span for TF motif enrichment
write.anchorBed.from.id <- function(x){

  df <- get(x)
  prefix <- str_replace(x,".domi.SuperAnchor","")%>%
            str_replace_all(.,"\\.","_")
  dir <- "../03_Super_Anchor/04_Dominate_superAnchor_bed/"

 df%>%rowwise%>%
    mutate(V1=str_split(anchor.id,"_",simplify = T)[1],
           V2=str_split(anchor.id,"_",simplify = T)[2]%>%as.numeric(),
           V3=str_split(anchor.id,"_",simplify = T)[3]%>%as.numeric()
           )%>%
    select(V1,V2,V3,anchor.id)%>%
    mutate(V2=V2-500,V3=V3+500)%>%
    arrange(V1,V2,V3)%>%
    write.table(file = paste0(dir,prefix,"_dominate_superAnchor.bed"),quote = F,col.names = F,row.names = F,sep = "\t")
}

ls(pattern = "domi.SuperAnchor")%>%lapply(.,write.anchorBed.from.id )

##Add the classifications for super anchors and dominate super anchors to the tbl object.
##And export the tbl objects to Rdata so that they can be load for other analysis.
lapply(c("D","N","RE"), function(x){

  df <- get(paste0(x,".tbl"))

 df %>% mutate(is.superAnchor.NC=case_when(anchor.id %in% N.superAnchors$anchor.id ~ "Yes",TRUE ~ "NO"),
                         is.superAnchor.DS=case_when(anchor.id %in% D.superAnchors$anchor.id ~ "Yes",TRUE ~ "NO"),
                         is.superAnchor.RE=case_when(anchor.id %in% RE.superAnchors$anchor.id ~ "Yes",TRUE ~ "NO"),
                         is.domi.superAnchor.NC=case_when(anchor.id %in% N.domi.SuperAnchor$anchor.id ~ "Yes",TRUE ~ "NO"),
                         is.domi.superAnchor.DS=case_when(anchor.id %in% D.domi.SuperAnchor$anchor.id ~ "Yes",TRUE ~ "NO"),
                         is.domi.superAnchor.RE=case_when(anchor.id %in% RE.domi.SuperAnchor$anchor.id ~ "Yes",TRUE ~ "NO"))%>%
   assign(paste0(x,".tbl"),.,envir = .GlobalEnv)

})

save(D.tbl, file = "DS.tbl.Rdata")
save(N.tbl, file = "NC.tbl.Rdata")
save(RE.tbl, file = "RE.tbl.Rdata")

##Add info on Hub-anchors for supplementary data
##Export the tbl to data frame at this step
lapply(c("D","N","RE"), function(x){

  df <- get(paste0(x,".tbl"))

  df %<>% mutate(is.hubAnchor.NC=case_when(anchor.id %in% N.HubAnchor$anchor.id ~ "Yes",TRUE ~ "NO"),
                is.hubAnchor.DS=case_when(anchor.id %in% D.HubAnchor$anchor.id ~ "Yes",TRUE ~ "NO"),
                is.hubAnchor.RE=case_when(anchor.id %in% RE.HubAnchor$anchor.id ~ "Yes",TRUE ~ "NO"))
  #node attribute
  df %N>% as.data.frame()%>%
    write.xlsx(paste0("../01_PPI_network_analysis/01_PPI_network_files/",x,"_PPInetwork.node.xlsx"))

  #edge attribute
  df %E>%as.data.frame()%>%
    write.xlsx(paste0("../01_PPI_network_analysis/01_PPI_network_files/",x,"_PPInetwork.edge.xlsx"))
})


#########===============================================================================================#########
#########Analyze the superiority of the condition-dominate super-anchors.
#########===============================================================================================#########
###Assess if the adjacent genes of the dominate super-anchors have superiority in expression or FC compared to other
###Hub-anchors or normal anchors.
column.in <- c("anchor.id","RW.group","DS.TPM.anchor","NC.TPM.anchor","RE.TPM.anchor","log2FC_DS_vs_NC.anchor",
               "log2FC_RE_vs_DS.anchor","avgTPM.DS.RWgroup","avgTPM.NC.RWgroup","avgTPM.RE.RWgroup")

graph.D.TPM.domi_vs_hub <- rbind(
D.tbl%>%as.data.frame%>%filter(is.domi.superAnchor.DS=="Yes")%>%
  select(all_of(column.in))%>%mutate(cat="domi.SPE"),
D.tbl%>%as.data.frame%>%filter(anchor.id %in% D.HubAnchor$anchor.id,is.domi.superAnchor.DS=="NO")%>%
  select(all_of(column.in))%>%mutate(cat="hub")
)%>%ggplot(aes(y=log2(DS.TPM.anchor),x=cat))+
    geom_boxplot()+
    ggsignif::geom_signif(comparisons = list(c("domi.SPE","hub")))

graph.N.TPM.domi_vs_hub <- rbind(
  N.tbl%>%as.data.frame%>%filter(is.domi.superAnchor.DS=="Yes")%>%
    select(all_of(column.in))%>%mutate(cat="domi.SPE"),
  N.tbl%>%as.data.frame%>%filter(anchor.id %in% D.HubAnchor$anchor.id,is.domi.superAnchor.DS=="NO")%>%
    select(all_of(column.in))%>%mutate(cat="hub")
)%>%ggplot(aes(y=log2(NC.TPM.anchor),x=cat))+
  geom_boxplot()+
  ggsignif::geom_signif(comparisons = list(c("domi.SPE","hub")))

graph.RE.TPM.domi_vs_hub <- rbind(
  D.tbl%>%as.data.frame%>%filter(is.domi.superAnchor.DS=="Yes")%>%
    select(all_of(column.in))%>%mutate(cat="domi.SPE"),
  D.tbl%>%as.data.frame%>%filter(anchor.id %in% D.HubAnchor$anchor.id,is.domi.superAnchor.DS=="NO")%>%
    select(all_of(column.in))%>%mutate(cat="hub")
)%>%ggplot(aes(y=log2(RE.TPM.anchor),x=cat))+
  geom_boxplot()+
  ggsignif::geom_signif(comparisons = list(c("domi.SPE","hub")))

pdf("../01_PPI_network_analysis/02_PPI_network_profiles/08_direct_PPI_TPM_domi_vs_hub.pdf",height = 4,width = 6,bg=bg)
ggpubr::ggarrange(graph.D.TPM.domi_vs_hub,graph.N.TPM.domi_vs_hub,graph.RE.TPM.domi_vs_hub,nrow = 1,
                  common.legend = T,align = "hv")
dev.off()
#The SPEs have no superiority on the activation of nearby genes.

#########===============================================================================================#########
####Assess if the dominate SPE can integrate more contacts by the dominate loops than ordinary hub anchors.
#DS
RWgroup.by.DS.dominate.SPE <- D.tbl%>%filter(anchor.id %in% D.domi.SuperAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

RWgroup.by.DS.hub <- D.tbl%>%filter(!(anchor.id %in% D.domi.SuperAnchor$anchor.id),anchor.id %in% D.HubAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

gg.DS.domiPPI.domiSPEvsHub <-
rbind(
#by dominate SPEs.
D.tbl %>% filter(RW.group %in% RWgroup.by.DS.dominate.SPE)%>%
          morph(to_split,RW.group,split_by = "nodes")%N>%
          mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
          filter(loop.domination == "DS_dominate")%N>%
          mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
          unmorph()%>%
          mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
          select(RW.group,dominate.PPI.content)%>%
          as.data.frame()%>%
          mutate(condition="DS",anchor.type="Dominate SPE")%>%
          unique(),
#by hubs.
D.tbl %>% filter(RW.group %in% RWgroup.by.DS.hub,!(RW.group %in% RWgroup.by.DS.dominate.SPE))%>%
          morph(to_split,RW.group,split_by = "nodes")%N>%
          mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
          filter(loop.domination == "DS_dominate")%N>%
          mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
          unmorph()%>%
          mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
          select(RW.group,dominate.PPI.content)%>%
          as.data.frame()%>%
          mutate(condition="DS",anchor.type="hub")%>%
          unique()
)%>%
  #filter(Dominate.PPI.in.RWgroup != max(Dominate.PPI.in.RWgroup))%>% #remove the largest group, which is a outlier.
ggplot(aes(x=anchor.type,y=dominate.PPI.content))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))

#NC
RWgroup.by.NC.dominate.SPE <- N.tbl%>%filter(anchor.id %in% N.domi.SuperAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

RWgroup.by.NC.hub <- N.tbl%>%filter(!(anchor.id %in% N.domi.SuperAnchor$anchor.id),anchor.id %in% N.HubAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

gg.NC.domiPPI.domiSPEvsHub <-
  rbind(
    #by dominate SPEs.
    N.tbl %>% filter(RW.group %in% RWgroup.by.NC.dominate.SPE)%>%
      morph(to_split,RW.group,split_by = "nodes")%N>%
      mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
      filter(loop.domination == "NC_dominate")%N>%
      mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
      unmorph()%>%
      mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
      #filter(PPI.in.RWgroup>2)%>%
      select(RW.group,dominate.PPI.content)%>%
      as.data.frame()%>%
      mutate(condition="NC",anchor.type="Dominate SPE")%>%
      unique(),
    #by hubs.
    N.tbl %>% filter(RW.group %in% RWgroup.by.NC.hub,!(RW.group %in% RWgroup.by.NC.dominate.SPE))%>%
      morph(to_split,RW.group,split_by = "nodes")%N>%
      mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
      filter(loop.domination == "NC_dominate")%N>%
      mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
      unmorph()%>%
      mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
      #filter(PPI.in.RWgroup>6)%>%
      select(RW.group,dominate.PPI.content)%>%
      as.data.frame()%>%
      mutate(condition="NC",anchor.type="hub")%>%
      unique()
  )%>%
  #filter(dominate.PPI.content != max(dominate.PPI.content))%>% #remove the largest group, which is a outlier.
  ggplot(aes(x=anchor.type,y=dominate.PPI.content))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))

#RE
RWgroup.by.RE.dominate.SPE <- RE.tbl%>%filter(anchor.id %in% RE.domi.SuperAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

RWgroup.by.RE.hub <- RE.tbl%>%filter(!(anchor.id %in% RE.domi.SuperAnchor$anchor.id),anchor.id %in% RE.HubAnchor$anchor.id)%>%
  as.data.frame()%>%pull(RW.group)

gg.RE.domiPPI.domiSPEvsHub <-
  rbind(
    #by dominate SPEs.
    RE.tbl %>% filter(RW.group %in% RWgroup.by.RE.dominate.SPE)%>%
      morph(to_split,RW.group,split_by = "nodes")%N>%
      mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
      filter(loop.domination == "RE_dominate")%N>%
      mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
      unmorph()%>%
      mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
      select(RW.group,dominate.PPI.content)%>%
      as.data.frame()%>%
      mutate(condition="RE",anchor.type="Dominate SPE")%>%
      unique(),
    #by hubs.
    RE.tbl %>% filter(RW.group %in% RWgroup.by.RE.hub,!(RW.group %in% RWgroup.by.RE.dominate.SPE))%>%
      morph(to_split,RW.group,split_by = "nodes")%N>%
      mutate(PPI.in.RWgroup=dim(.E())[1])%E>%
      filter(loop.domination == "RE_dominate")%N>%
      mutate(Dominate.PPI.in.RWgroup=dim(.E())[1])%>%
      unmorph()%>%
      mutate(dominate.PPI.content=Dominate.PPI.in.RWgroup/PPI.in.RWgroup)%>%
      select(RW.group,dominate.PPI.content)%>%
      as.data.frame()%>%
      mutate(condition="RE",anchor.type="hub")%>%
      unique()
  )%>%
  #filter(Dominate.PPI.in.RWgroup != max(Dominate.PPI.in.RWgroup))%>% #remove the largest group, which is a outlier.
  ggplot(aes(x=anchor.type,y=dominate.PPI.content))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))


ggpubr::ggarrange(gg.DS.domiPPI.domiSPEvsHub,gg.NC.domiPPI.domiSPEvsHub,gg.RE.domiPPI.domiSPEvsHub,
                  align = "hv",ncol = 1,common.legend = T)

###View genes in the RW groups by dominate SPE
###View genes in the RW groups by dominate SPE
D.tbl%N>%filter(RW.group %in% RWgroup.by.DS.dominate.SPE)%>%
  as.data.frame()%>%select(name,RW.group,degree.dens)%>%
  unique()%>%filter(name != "")%>%
  left_join(D.tbl%N>%filter(anchor.id %in% D.domi.SuperAnchor$anchor.id)%>%as.data.frame()%>%select(anchor.id,RW.group))%>%na.omit

N.tbl%N>%filter(RW.group %in% RWgroup.by.NC.dominate.SPE)%>%
  as.data.frame()%>%select(name,RW.group,degree.dens)%>%
  unique()%>%filter(name != "")%>%
  left_join(D.tbl%N>%filter(anchor.id %in% N.domi.SuperAnchor$anchor.id)%>%as.data.frame()%>%select(anchor.id,RW.group))%>%na.omit

#########===============================================================================================#########
####Assess if the RW group integrated by the dominate SPE have superiority in gene expression level or changes compared
####to the other groups lead by the normal hub genes.
##NC
gg.NC.FC.domiSPEvsHub <-
  rbind(
    #by dominate SPEs.
    N.tbl %>% filter(RW.group %in% RWgroup.by.NC.dominate.SPE)%E>%
      filter(loop.domination == "NC_dominate")%E>%
      mutate(log2FC_DS_vs_NC=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_DS_vs_NC)%>%
      mutate(condition="NC",anchor.type="Dominate SPE")%>%
      unique(),
    #by hubs.
    N.tbl %>% filter(RW.group %in% RWgroup.by.NC.hub)%E>%
      filter(loop.domination == "NC_dominate")%E>%
      mutate(log2FC_DS_vs_NC=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_DS_vs_NC)%>%
      mutate(condition="NC",anchor.type="hub")%>%
      unique()
  )%>%
  ggplot(aes(x=anchor.type,y=-log2FC_DS_vs_NC))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  ylab("NC-Dominate PPI genes log2FC_NC_vs_DS")+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))


##DS
gg.DS.FC.domiSPEvsHub <-
  rbind(
    #by dominate SPEs.
    D.tbl %>% filter(RW.group %in% RWgroup.by.DS.dominate.SPE)%E>%
      filter(loop.domination == "DS_dominate")%E>%
      mutate(log2FC_DS_vs_NC=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_DS_vs_NC)%>%
      mutate(condition="DS",anchor.type="Dominate SPE")%>%
      unique(),
    #by hubs.
    D.tbl %>% filter(RW.group %in% RWgroup.by.DS.hub)%E>%
      filter(loop.domination == "DS_dominate")%E>%
      mutate(log2FC_DS_vs_NC=(log2FC_DS_vs_NC.anchor.1+log2FC_DS_vs_NC.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_DS_vs_NC)%>%
      mutate(condition="DS",anchor.type="hub")%>%
      unique()
  )%>%
  ggplot(aes(x=anchor.type,y=log2FC_DS_vs_NC))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  ylab("DS-Dominate PPI genes log2FC_DS_vs_NC")+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))

##DS
gg.RE.FC.domiSPEvsHub <-
  rbind(
    #by dominate SPEs.
    RE.tbl %>% filter(RW.group %in% RWgroup.by.RE.dominate.SPE)%E>%
      filter(loop.domination == "RE_dominate")%E>%
      mutate(log2FC_RE_vs_DS=(log2FC_RE_vs_DS.anchor.1+log2FC_RE_vs_DS.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_RE_vs_DS)%>%
      mutate(condition="RE",anchor.type="Dominate SPE")%>%
      unique(),
    #by hubs.
    RE.tbl %>% filter(RW.group %in% RWgroup.by.RE.hub)%E>%
      filter(loop.domination == "RE_dominate")%E>%
      mutate(log2FC_RE_vs_DS=(log2FC_RE_vs_DS.anchor.1+log2FC_RE_vs_DS.anchor.2)/2)%E>%
      as.data.frame()%>%
      select(log2FC_RE_vs_DS)%>%
      mutate(condition="RE",anchor.type="hub")%>%
      unique()
  )%>%
  ggplot(aes(x=anchor.type,y=log2FC_RE_vs_DS))+
  geom_violin(trim = F,draw_quantiles = F)+
  geom_boxplot(outlier.shape = NA,width=0.3)+
  ylab("RE-Dominate PPI genes log2FC_RE_vs_DS")+
  geom_signif(comparisons = list(c("Dominate SPE","hub")))


ggpubr::ggarrange(gg.NC.FC.domiSPEvsHub,gg.DS.FC.domiSPEvsHub,gg.RE.FC.domiSPEvsHub,
                  align = "hv",ncol = 1,common.legend = T)


###Plot the dominate PPI content and dominate PPI log2FC together.
pdf("../01_PPI_network_analysis/02_PPI_network_profiles/09_domiPPI_FC_RWgroup_domiSPE_vs_hub.pdf",height = 12,width = 3,bg=bg)
ggpubr::ggarrange(gg.DS.domiPPI.domiSPEvsHub,gg.NC.domiPPI.domiSPEvsHub,gg.RE.domiPPI.domiSPEvsHub,
                  gg.NC.FC.domiSPEvsHub,gg.DS.FC.domiSPEvsHub,gg.RE.FC.domiSPEvsHub,
                  align = "hv",ncol = 1,common.legend = T)
dev.off()

#########===============================================================================================#########
#########Plot PPI network of each condition.
#########===============================================================================================#########
#We can see that the anchor network is composed of a large connected subnetwork (involoving thousands of anchors)
#and numerous isolated tiny interacting pairs with only 2-10 anchors.
#We plot the largest subnetwork by centrality layout to highlight the connecting hub anchor and the hierarchical distribution of
#anchor interacting frequency.

node.column <- c("anchor.id","degree","degree.dens","represent.gene","DS.TPM.anchor","NC.TPM.anchor",
              "RE.TPM.anchor","log2FC_DS_vs_NC.anchor","log2FC_RE_vs_DS.anchor",
              "log2FC_RE_vs_NC.anchor","mfuzz.cluster.anchor","membership.mfuzz.1","membership.mfuzz.2","membership.mfuzz.3",
              "membership.mfuzz.4","membership.mfuzz.5","membership.mfuzz.6","MSU7_ID","name","is.superAnchor.NC",
              "is.superAnchor.DS","is.superAnchor.RE","is.domi.superAnchor.NC","is.domi.superAnchor.DS","is.domi.superAnchor.RE")

tbl <- graph_join(D.tbl%N>%select(all_of(node.column))%>%rename(degree.DS=degree,degree.dens.DS=degree.dens),
                  N.tbl%N>%select(all_of(node.column))%>%rename(degree.NC=degree,degree.dens.NC=degree.dens)
                  )%>%
        graph_join(RE.tbl%N>%select(all_of(node.column))%>%rename(degree.RE=degree,degree.dens.RE=degree.dens))

#DS.anchor.to.domi.PPI <- D.tbl%E>%filter(loop.domination=="DS_dominate")%N>%filter(!node_is_isolated())%N>%pull(anchor.id)
#NS.anchor.to.domi.PPI <- N.tbl%E>%filter(loop.domination=="NC_dominate")%N>%filter(!node_is_isolated())%N>%pull(anchor.id)
#RE.anchor.to.domi.PPI <- RE.tbl%E>%filter(loop.domination=="RE_dominate")%N>%filter(!node_is_isolated())%N>%pull(anchor.id)

#tbl %<>% mutate(is.DS.domi.superAnchor = ifelse(anchor.id %in% D.superAnchors$anchor.id,"Yes","NO"),
#                is.NC.domi.superAnchor = ifelse(anchor.id %in% N.superAnchors$anchor.id,"Yes","NO"),
#                is.RE.domi.superAnchor = ifelse(anchor.id %in% RE.superAnchors$anchor.id,"Yes","NO"))

#Set the anchor genes to be labeled.
#filter super anchor nodes
tbl%>%as.data.frame()%>%arrange(desc(log2FC_DS_vs_NC.anchor))%>%filter(name != "")%>%
  filter(anchor.id %in% D.superAnchors$anchor.id)%>%select(name,degree.DS,degree.NC,degree.RE)%>%unique()

#filter node connected by dominate loops.
tbl%E>%filter(loop.domination=="DS_dominate")%N>%filter(!node_is_isolated())%>%
  as.data.frame()%>%
  arrange(desc(log2FC_DS_vs_NC.anchor))%>%filter(name != "")%>%
  #filter(anchor.id %in% D.superAnchors$anchor.id)%>%
  select(name,log2FC_DS_vs_NC.anchor)%>%unique()

tbl%E>%filter(loop.domination=="NC_dominate")%N>%filter(!node_is_isolated())%>%
  as.data.frame()%>%
  arrange(log2FC_DS_vs_NC.anchor)%>%filter(name != "")%>%
  #filter(anchor.id %in% D.superAnchors$anchor.id)%>%
  select(name,log2FC_DS_vs_NC.anchor)%>%unique()


#to label super-anchors
anchor.label <- c("DST","OsjDHN6|RAB16C","OsjDHN7|RAB16B","OsjDHN5","LEA16","OsjDHN6|RAB16C",
                  "HSFB2C")

#to label DS dominate PPI loop-associated genes.
#"OsHSF6|OsHsfA2c|OsHsfA6a" site is an example for intra-CID interaction reprogramming.
anchor.label2 <- c("DST","OsPYL6",                #these are DS-suppressed genes
                   "OsjDHN6|RAB16C","OsjDHN7|RAB16B","OsjDHN5","LEA16","OsjDHN6|RAB16C",
                  "HSFB2C","OsABIL3|OsPP2C50","DREB1C","OSHB4|OsHox32","OsHSF6|OsHsfA2c|OsHsfA6a",
                  "OsHSF6|OsHsfA2c|OsHsfA6a","HSFA1","OsJAZ9","SNAC1|OsNAC19|OsNAC9","OsMKK10-2|MPKK10.2")

#Plot networks
#The GW network is too large for clear view.
#We plot the networks within the selected chromosomes.
#Calculate the layout for total network. This will make the subsequent plotting easier.
tbl.chr01 <- tbl%>%filter(str_split(anchor.id,"_",simplify = T)[,1]=="Chr01")
tbl.chr03 <- tbl%>%filter(str_split(anchor.id,"_",simplify = T)[,1]=="Chr03")
tbl.chr11 <- tbl%>%filter(str_split(anchor.id,"_",simplify = T)[,1]=="Chr11")
tbl.chr04 <- tbl%>%filter(str_split(anchor.id,"_",simplify = T)[,1]=="Chr04")

#tbl.weights <- tbl%E>%as.data.frame%>%select(ipet.counts)%>%pull()
#tbl.xy <- layout_with_fr(tbl,weights=tbl.weights,grid="nogrid")
set.seed(100)
tbl.chr01.weights <-tbl.chr01%E>%as.data.frame%>%select(ipet.counts)%>%pull()
tbl.chr01.xy <- layout_with_fr(tbl.chr01,grid="nogrid")

tbl.chr03.weights <-tbl.chr03%E>%as.data.frame%>%select(ipet.counts)%>%pull()
tbl.chr03.xy <- layout_with_fr(tbl.chr03,grid="nogrid")

tbl.chr11.weights <-tbl.chr11%E>%as.data.frame%>%select(ipet.counts)%>%pull()
tbl.chr11.xy <- layout_with_fr(tbl.chr11,grid="nogrid")

tbl.chr04.weights <-tbl.chr04%E>%as.data.frame%>%select(ipet.counts)%>%pull()
tbl.chr04.xy <- layout_with_fr(tbl.chr04,grid="nogrid")


#The layout are used for the subsequent plot, and node positions will be fixed among conditions.

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
  geom_node_point(aes(filter=is.superAnchor.DS=="Yes"),
                  shape=23,size=8,fill=col.DS)+
  geom_node_point(aes(filter=is.domi.superAnchor.DS=="Yes"),
                  shape=22,size=8,fill=col.DS)+
  scale_fill_manual(values=c("1"="lightgreen",
                             "2"="pink",
                             "3"="pink",
                             "4"="cyan",
                             "5"="lightgreen",
                             "6"="skyblue",
                             "NA"="grey50"))+
  geom_node_text(aes(label=name,filter=name %in% anchor.label2),size=5)+
  theme_graph(base_family= "Arial")+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(title="")
dev.off()

#--plot NC
pdf("../01_PPI_network_analysis/03_PPI_network_plot/NC.allPPIs.Chr11.pdf",height = 15,width = 15,bg=bg)
ggraph(tbl.chr11,layout = "manual",x=tbl.chr11.xy[,1],y=tbl.chr11.xy[,2])+
  geom_edge_density(aes(fill=loop.domination,
                        filter=loop.in.NC==1))+
  scale_edge_fill_manual(values=c("other"=bg,
                                  "Conserved"=bg,
                                  "RE_dominate"=col.RE,
                                  "NC_dominate"=col.NC,
                                  "DS_dominate"=col.DS))+
  geom_edge_link(aes(edge_color=loop.domination,
                     edge_width=loop.domination,
                     filter=loop.in.NC==1),
                 check_overlap =T)+
  scale_edge_width_manual(values=c("other"=0.8,
                                   "Conserved"=0.8,
                                   "RE_dominate"=2,
                                   "NC_dominate"=2,
                                   "DS_dominate"=2))+
  scale_edge_color_manual(values=c("other"="grey75",
                                   "Conserved"="grey70",
                                   "RE_dominate"=bg,
                                   "NC_dominate"=col.NC,
                                   "DS_dominate"=bg))+
  geom_node_point(aes(fill=as.factor(mfuzz.cluster.anchor),
                      size=membership.mfuzz.1),
                  shape=21)+
  scale_size_continuous(range = c(0.1,8))+
  #highlight DS-dominate super-anchors
  geom_node_point(aes(filter=is.superAnchor.NC=="Yes"),
                  shape=23,size=8,fill=col.NC)+
  geom_node_point(aes(filter=is.domi.superAnchor.NC=="Yes"),
                  shape=22,size=8,fill=col.NC)+
  scale_fill_manual(values=c("1"="lightgreen",
                             "2"="pink",
                             "3"="pink",
                             "4"="cyan",
                             "5"="lightgreen",
                             "6"="skyblue",
                             "NA"="grey50"))+
  geom_node_text(aes(label=name,filter=name %in% anchor.label2),size=5)+
  theme_graph(base_family= "Arial")+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(title="")
dev.off()

#--plot RE
pdf("../01_PPI_network_analysis/03_PPI_network_plot/RE.allPPIs.Chr04.pdf",height = 20,width = 20,bg=bg)
ggraph(tbl.chr04,layout = "manual",x=tbl.chr04.xy[,1],y=tbl.chr04.xy[,2])+
  geom_edge_density(aes(fill=loop.domination,
                        filter=loop.in.RE==1))+
  scale_edge_fill_manual(values=c("other"=bg,
                                  "Conserved"=bg,
                                  "RE_dominate"=col.RE,
                                  "NC_dominate"=col.NC,
                                  "DS_dominate"=col.DS))+
  geom_edge_link(aes(edge_color=loop.domination,
                     edge_width=loop.domination,
                     filter=loop.in.RE==1),
                 check_overlap =T)+
  scale_edge_width_manual(values=c("other"=0.8,
                                   "Conserved"=0.8,
                                   "RE_dominate"=2,
                                   "NC_dominate"=2,
                                   "DS_dominate"=2))+
  scale_edge_color_manual(values=c("other"="grey75",
                                   "Conserved"="grey70",
                                   "RE_dominate"=col.RE,
                                   "NC_dominate"=bg,
                                   "DS_dominate"=bg))+
  geom_node_point(aes(fill=as.factor(mfuzz.cluster.anchor),
                      size=membership.mfuzz.1),
                  shape=21)+
  scale_size_continuous(range = c(0.1,8))+
  #highlight DS-dominate super-anchors
  geom_node_point(aes(filter=is.superAnchor.RE=="Yes"),
                  shape=23,size=8,fill=col.RE)+
  geom_node_point(aes(filter=is.domi.superAnchor.RE=="Yes"),
                  shape=22,size=8,fill=col.RE)+
  scale_fill_manual(values=c("1"="lightgreen",
                             "2"="pink",
                             "3"="pink",
                             "4"="cyan",
                             "5"="lightgreen",
                             "6"="skyblue",
                             "NA"="grey50"))+
  geom_node_text(aes(label=name,filter=name %in% anchor.label2),size=5)+
  theme_graph(base_family= "Arial")+
  theme(legend.position = "none")+
  coord_fixed()+
  labs(title="")
dev.off()
