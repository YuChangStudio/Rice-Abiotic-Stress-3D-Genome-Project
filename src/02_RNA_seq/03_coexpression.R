#######==================================================================================================#######
#######In this section, we characterize co-expression genes by three methods:
#######1. Dynamic time warp;
#######2. mfuzz;
#######3. Pearson correlation coefficiency;
#######4. GENIE3 algorithm (script derived from: https://doi.org/10.1038/s41467-021-26165-3).
#######==================================================================================================#######

#######==================================================================================================#######
#######Load configs
#######==================================================================================================#######

library(dtwclust)
library(ggplot2)
library(ggbeeswarm)
library(openxlsx)
library(stringr)
library(dplyr)
library(magrittr)

theme(
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

bg="transparent"

#GENIE3 edge weight calculation script.
#source("RS.Get.Weight.Matrix.CY.R")

#Use TPM values of all expressed genes for the analysis
TPM.exp <- read.delim("../01_in_matrix/gene.tpm",header = T)%>%
  filter(MH_D_1+MH_D_2 > 2 | MH_N_1+MH_N_2 > 2 | MH_RE_1+MH_RE_2 > 2)%>%
  mutate(DS=(MH_D_1+MH_D_2)/2,NC=(MH_N_1+MH_N_2)/2,RE=(MH_RE_1+MH_RE_2)/2)%>%
  select(gene_id,NC,DS,RE)

rownames(TPM.exp) <- TPM.exp$gene_id

TPM.exp <- TPM.exp%>%select(-gene_id)


ricefungene <- read.xlsx("~/genome_anno/MH63RS3/MH63RS3_encode_Sep2022.xlsx")

#######==================================================================================================#######
#######Identify DS-RE treatments drived co-expression groups by mfuzz
#######==================================================================================================#######
library(Mfuzz)
library(viridis)
#Standardize the TPM data
in.data <- new('ExpressionSet',exprs = TPM.exp%>%as.matrix)

in.data <- standardise(in.data)

#Optimize the m value to avoid random clustering
m <- mestimate(in.data)

#we set 6 clusters and perform clustering
set.seed(100)
fuzz.res <- mfuzz(in.data, c = 6, m = m)

fuzz.res$cluster%>%table()

#Assess membership between gene and cluster.
head(fuzz.res$membership)

#Assess the expressing trend of each cluster
#Plot all genes and core genes (mem >0.75) respectively.
#Depth of the line color represent membership values.
pdf("../02_mfuzz/CoExpression.trend.min.mem=0.pdf",height = 6,width = 8,bg="transparent")
mfuzz.plot(in.data, cl = fuzz.res, mfrow = c(2, 3),
           time.labels = c("NC","DS","RE"),
           new.window=F,
           min.mem=0,
           colo=rev(inferno(60)))
dev.off()

pdf("../02_mfuzz/CoExpression.trend.min.mem=0.5.pdf",height = 6,width = 8,bg="transparent")
mfuzz.plot(in.data, cl = fuzz.res, mfrow = c(2, 3),
           time.labels = c("NC","DS","RE"),
           new.window=F,
           min.mem=0.5,
           colo=rev(inferno(60)))
dev.off()

#Output clustering info
fuzz.res$cluster%>%as.data.frame()%>%rename(mfuzz.cluster=".")%>%
  mutate(gid=rownames(.))%>%
  merge(fuzz.res$membership,by="row.names",all.x=T)%>%
  left_join(.,ricefungene%>%select(MH63RS3_ID,MSU7_ID,name,Title),by=c("gid"="MH63RS3_ID"))%>%
  select(-Row.names)%>%
  write.xlsx(file = "../02_mfuzz/ExpressedGenes.mfuzzModuleMemship.xlsx",x=.)

#######==================================================================================================#######
########Investigate the relationships between each module.
#######==================================================================================================#######

pdf("../02_mfuzz/mfuzz.module.overlap015.PCA.pdf",height = 6,width = 6,bg="transparent")
overlap.plot(fuzz.res,overlap(fuzz.res),thres = 0.15)
dev.off()

pdf("../02_mfuzz/mfuzz.module.overlap02.PCA.pdf",height = 6,width = 6,bg="transparent")
overlap.plot(fuzz.res,overlap(fuzz.res),thres = 0.2)
dev.off()

pdf("../02_mfuzz/mfuzz.module.overlap01.PCA.pdf",height = 6,width = 6,bg="transparent")
overlap.plot(fuzz.res,overlap(fuzz.res),thres = 0.1)
dev.off()

#######==================================================================================================#######
#######Integrate the co-expression result with condition-dominate PPI relations
#######==================================================================================================#######

detach("package:dplyr")
library(multcomp)
library(dplyr)

#Load the representative gene info for DS dominate PPI.
#Set 1, 2, 3, 6 refer to NC-dominate, consensous, RE-dominate, and DS-dominate, respectively.
PPI.intersects <- read.xlsx("../../../22_diffloop_analysis/01_diffLoops_treatments/04_diffPPI_list/PPI.intersectionSet.NC.DS.RE.xlsx")

#Get the PPI genes of each set
lapply(seq(1:7), function(x){

  setx <- paste0("set",x)

  PPIgene.set <- PPI.intersects%>%filter(set==paste0(setx))

  PPIgene.set <- c(PPIgene.set$represent.gene.1,PPIgene.set$represent.gene.2)%>%unique()

  assign(paste0("PPIgene.set",x),PPIgene.set,envir = .GlobalEnv)

  #Get membership values for each set

  mem <- fuzz.res$membership%>%
         as.data.frame()%>%
         mutate(gid=rownames(.))%>%
         mutate(cat=case_when(gid%in%PPIgene.set ~ "in module",
                              TRUE ~ "other"),
                set=paste0(setx))%>%
          reshape2::melt()%>%
          rename(module=variable,membership=value)

  assign(paste0("mem.set",x),mem,envir = .GlobalEnv)

})


mem.all.sets <- rbind(mem.set1,mem.set2,mem.set3,mem.set4,mem.set5,mem.set6,mem.set7)

#Statistics on the mfuzz co-expression memberships of the DS-dominate PPI genes.

pdf("../03_coExp_wi_PPI/Dominate_PPI_coExp_memship.box.pdf",height = 8, width = 8, bg="transparent")
ggplot(data = mem.all.sets, aes(x=module,y=log2(membership),fill=cat))+
geom_boxplot(position='dodge', outlier.shape = NA)+
facet_wrap(.~set,ncol=3,nrow = 3,scales = "free")
#geom_violin(trim = F)
#geom_half_violin(draw_quantiles = F,trim = T)
dev.off()

#Anova with HSD
lapply(seq(1:7), function(x){

  mem.hsd<- mem.all.sets%>%
    filter(set == paste0("set",x))%>%
    mutate(col=str_c(module,cat,sep = "_"))%>%
    mutate(col=as.factor(col))%>%
    mutate(log2membership=log2(membership))%>%
    lm(log2membership ~ col,data = .)%>%glht(.,linfct = mcp(col = "Tukey"))%>%cld(decreasing=T)

  assign(paste0("mem.hsd.set",x),mem.hsd,envir = .GlobalEnv)
})


#label the significances in AI
write.table(membership.PPI.hsd$mcletters$Letters
            ,"../03_coExp_wi_PPI/membership.DsDomPPI.hsd.txt",quote = F)


#######==================================================================================================#######
#######Plot the GO results for each co-expression module
#######==================================================================================================#######
#Load GO results by agriGO V2
lapply(seq(1:6), function(x){

  go <- read.delim(paste0("../02_mfuzz/GO/02_selected/module.",x,".GO.txt"))%>%
        select(GO_acc,term_type,Term,queryitem,querytotal,bgitem,bgtotal,FDR)%>%
        filter(FDR < 0.01)%>%
        mutate(module=paste0("module",x))%>%
        mutate(term_type = factor(term_type,level=c("P","F","C"),ordered = T),
               Term = factor(Term,level=rev(Term),ordered = T))

  assign(paste0("GO.module.",x),go,envir = .GlobalEnv)

})

GO.all <- rbind(GO.module.1,GO.module.2,GO.module.3,GO.module.4,GO.module.5,GO.module.6)


#Plot GO result
library(viridis)
library(ggpubr)


pm1 <- ggplot(GO.module.1) +
          geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
          facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
          scale_fill_viridis()

pm2 <- ggplot(GO.module.2) +
          geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
  facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
  scale_fill_viridis()


pm3 <- ggplot(GO.module.3) +
        geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
  facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
  scale_fill_viridis()

pm4 <- ggplot(GO.module.4) +
        geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
  facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
  scale_fill_viridis()


pm5 <- ggplot(GO.module.5) +
        geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
  facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
  scale_fill_viridis()


pm6 <- ggplot(GO.module.6) +
        geom_bar(aes(y=Term,x=queryitem,fill=FDR),stat = "identity")+
  facet_grid(rows = vars(term_type),scales = "free_y",space  = "free_y")+
  scale_fill_viridis()


pdf("../02_mfuzz/GO/03_GO_plot/GO.barplot.pdf",bg=bg,height = 8,width = 10 )
ggarrange(pm1,pm2,pm3,pm4,pm5,pm6,
          ncol = 2,nrow=3,
          align="hv",
          #labels = c("Module 1","Module 2","Module 3","Module 4","Module 5","Module 6"),
          common.legend = T)
dev.off()



#The co-expression membership to module 2 is significantly higher in DS-dominate PPI genes, while other
#modules do not exhibit such difference.

#Export the DS-dominate PPI genes in module 2:


#Export node and edge files for DS-dominate PPI network with modules annotated.
PPI.DSdom%>%rowwise()%>%
  mutate(chr.x=str_split(anchor1.id,pattern = "_",simplify = T)[1],
         chr.y=str_split(anchor2.id,pattern = "_",simplify = T)[1])%>%
  mutate(edge.chr=ifelse(chr.x==chr.y,chr.x,"inter_chr"))%>%
  select(flank_geneIds.x,flank_geneIds.y,loop.index,edge.chr)%>%
  write.table("../03_coExp_wi_PPI/DS_dominate_PPI.edge.txt",sep = "\t",quote = F,row.names = F)

rbind(PPI.DSdom%>%rename(node=flank_geneIds.x,module=mfuzz.cluster.x)%>%
        mutate(node.name=ifelse(name.x=="",node,name.x))%>%
        select(node,node.name,module),

      PPI.DSdom%>%rename(node=flank_geneIds.y,module=mfuzz.cluster.y)%>%
        mutate(node.name=ifelse(name.y=="",node,name.y))%>%
        select(node,node.name,module))%>%
  arrange(node)%>%
  write.table("../03_coExp_wi_PPI/DS_dominate_PPI.node.txt",sep = "\t",quote = F,row.names = F)


