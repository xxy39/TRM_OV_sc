
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(AUCell)
color.5group<-c("#bea0d2", "#e6194b" ,"#8eb166", "#66b3dd", "#4363d")

# read in seurat object
integerated<-readRDS("integrated.rds")

# markers for different T cell differentiation states, derived from  Anadon et al. 2022 Cancer Cell.
celltype.genelist<-list(Stemness=c("SELL","LEF1","CD28","CD27","CCR7","IL7R","CXCR5","TCF7","BACH2","JUNB","EGR1","KLF2"),
                        Effector=c("GZMK","GZMH","GZMB","PRF1","GNLY","IFNG","FASLG","FGFBP2"),
                        Proliferation=c("TOP2A","MKI67","CDK1","STMN1","DNMT1","MCM7"),
                        Exhaustion=c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","CD74","MIR155HG","CXCL13","LAYN","MYO7A","HLA-DRB1","HLA-DQA1",
                                     "TOX","TOX2","BATF","ETV1","ID2","ZNF683","RBPJ","TBX21","RUNX1","RUNX3","RUNX2","EZH2"))

# Transcription_Factor=c("TCF7","BACH2","EOMES","PRDM1","JUNB","EGR1","KLF2","TOX","TOX2","BATF","ETV1","ID2","ZNF683","RBPJ","TBX21","RUNX1","RUNX3","RUNX2","EZH2"))

# Gene set acvitity scores using AUCell
exprMatrix <- GetAssayData(object = integerated[['RNA']])
cells_rankings <- AUCell_buildRankings(as.matrix(exprMatrix), nCores=4, plotStats=TRUE)

# Calcualte AUCell scores for celltype.genelist. Cells with <1000 features were skipped
cells_AUC <- AUCell_calcAUC(celltype.genelist, cells_rankings,aucMaxRank=1000)

# assign scores to each cell
integerated$AUCell_Stemness<- getAUC(cells_AUC)[1,]
integerated$AUCell_Effector<- getAUC(cells_AUC)[2,]
integerated$AUCell_Proliferation<- getAUC(cells_AUC)[3,]
integerated$AUCell_Exhaustion<- getAUC(cells_AUC)[4,]

# generate dotplot for AUCell scores
dot<-DotPlot(integerated, features = c("AUCell_Stemness","AUCell_Effector","AUCell_Proliferation","AUCell_Exhaustion"))+ RotatedAxis()+scale_colour_gradient2(low ="#4575B4", mid = "#FFFFBF", high = "#D73027")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9))

print(dot)

# assign differentiation states to clusters based on AUC scores
# stemness: 0,2,4,7,15,20,21,22
# effector: 3,5,6,10,12,18,19
# Proliferative: 13, 23
# Exhaustive: 1,8,9,11,14,16,17

integerated$state<-rep("Stemness",ncol(integerated))
integerated$state[Idents(integerated)%in%c(3,5,6,10,12,18,19)]<-"Effector"
integerated$state[Idents(integerated)%in%c(13, 23)]<-"Proliferative"
integerated$state[Idents(integerated)%in%c(1,8,9,11,14,16,17)]<-"Exhaustive"

# calculate percentage of of TRMs in each cluster
integerated$Cluster<-Idents(integerated)
TRM<-integerated@meta.data %>% group_by(Cluster) %>% dplyr::summarize(TRM=length(which(CellType=="TRM"))/n(), Count=n())

# plot median stemness scores vs. median exhaustion scores for each cluster
# color the cluster by % of TRM cells
# size of the cluster corresponds to number of cells
AUC_median<-integerated@meta.data %>%
  group_by(Cluster) %>%
  dplyr::summarize(AUCell_Stemness_median = median(AUCell_Stemness, na.rm = TRUE),
                   AUCell_Exhaustion_median = median(AUCell_Exhaustion, na.rm = TRUE)) 

integerated.AUC.avg<- AUC_median %>% select(AUCell_Stemness_median,AUCell_Exhaustion_median ) %>% as.data.frame()
rownames(integerated.AUC.avg)<-levels(Idents(integerated))
integerated.AUC.avg.zcale<-as.data.frame(scale(integerated.AUC.avg))
integerated.AUC.avg.zcale$TRM.Pct<-TRM$TRM
integerated.AUC.avg.zcale$Count<-TRM$Count
integerated.AUC.avg.zcale$Cluster<-rownames(integerated.AUC.avg.zcale)

ggplot(integerated.AUC.avg.zcale, aes(x =AUCell_Exhaustion_median , y = AUCell_Stemness_median, label=Cluster)) + 
  geom_point(aes(fill = TRM.Pct,size = Count), pch=21, color="black", alpha=0.7)+
  scale_fill_gradient2(low ="#4575B4", mid = "#FFFFBF", high = "#D73027" , midpoint = 0.5, breaks=c(0.2,0.4,0.6,0.8),name = "TRM Percentage")+geom_text()+
  scale_radius(range = c(4,12),name = "Cell count")+ylab("Stemness Score")+xlab("Exhaustion Score")+
  theme_light()

