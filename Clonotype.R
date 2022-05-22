
library(ggplot2)
library(Seurat)
library(fossil)
library(reshape2)
library(readxl)    
color.5group<-c("#bea0d2", "#e6194b" ,"#8eb166", "#66b3dd", "#4363d8")


# load in integrated Seurat object
integerated<-readRDS("integrated.rds")

# subset integerated to 5 major groups

integerated.subset<-subset(integerated,cells=colnames(integerated)[integerated$Group%in%c("Stemlike_recir","Stemlike_TRM","Effector_TRM","Proliferative_TRM","Exhausted_TRM")])
# count number of cells in each clonotype and each group
clonotype.table<-integerated.subset@meta.data %>% select(clonotype_id.consensus,Group) %>% 
  tidyr::gather(clonotype_id.consensus,Group) %>% 
  dplyr::group_by(clonotype_id.consensus,Group) %>%  dplyr::summarise(n = n()) %>% 
  tidyr::spread(Group, n) %>% as.data.frame()
clonotype.table<-clonotype.table %>% filter(!is.na(clonotype_id.consensus))
rownames(clonotype.table)<-clonotype.table$clonotype_id.consensus
clonotype.table<-clonotype.table[,-1]
clonotype.table[is.na(clonotype.table)] <- 0

# calculate morisita.horn similarity scores of clonotypes between groups
groups.list<-levels(factor(integerated.subset$Group))

morisita.horn.matrix<-matrix(NA,ncol=length(groups.list), nrow=length(groups.list))

for(i in 1:(length(groups.list)-1)){
  morisita.horn.matrix[i,i]=1
  morisita.horn.matrix[i+1,i+1]=1
  
  for(j in (i+1):length(groups.list)){
    
    morisita.horn.matrix[i,j]=morisita.horn(clonotype.table[,groups.list[i]], clonotype.table[,groups.list[j]])
  }
}

colnames(morisita.horn.matrix)=groups.list
rownames(morisita.horn.matrix)=groups.list

# heatmap plot
melted_morisita.horn.matrix <- melt(morisita.horn.matrix, na.rm = TRUE)
melted_morisita.horn.matrix$Var1<-factor(melted_morisita.horn.matrix$Var1, levels=groups.list)
melted_morisita.horn.matrix$Var2<-factor(melted_morisita.horn.matrix$Var2, levels=groups.list)
morisita.horn.plot<-ggplot(data = melted_morisita.horn.matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.01, limit = c(0,1), space = "Lab", 
                       name="Morisita-Horn Similarity") +
  theme_light()+ 
  theme(axis.text.x = element_text( size = 12, angle = 45,hjust=1),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+ylab("Group")+xlab("Group")+
  scale_x_discrete(breaks=groups.list,labels=groups.list)+
  scale_y_discrete(breaks=groups.list,labels=groups.list)

morisita.horn.plot



# trajectory plot of individule clonotypes
TCR2plot<-readxl::read_excel("/Volumes/SSD_1TB/Dropbox/CICPT_1954/Paper_Figures/4tumors_scRNAseq/SelectedTCRs.xlsx",sheet = 1)
integerated.UMAP<-as.data.frame(Embeddings(object = integerated.subset, reduction = "umap"))
integerated.UMAP$Cluster<-integerated.subset$Cluster
integerated.UMAP$Group<-factor(integerated.subset$Group,levels=c("Stemlike_recir","Stemlike_TRM","Effector_TRM","Proliferative_TRM","Exhausted_TRM"))
integerated.UMAP$clonotype_id.consensus<-integerated.subset$clonotype_id.consensus

for(i in 1:nrow(TCR2plot)){
  integerated.UMAP$TCR<-ifelse(integerated.UMAP$clonotype_id.consensus==TCR2plot$clonotype_id.consensus[i],1,0)
  tcr.cell<-integerated.UMAP %>% filter(TCR==1 )
  p1<-ggplot(integerated.UMAP, aes(x=UMAP_1,y=UMAP_2))+geom_point(color="gray80",size=0.1)+
    geom_point(data=tcr.cell, aes(x=UMAP_1,y=UMAP_2, fill=Group),color="black", pch=21, size=1.5)+scale_fill_manual(name="Type",values=color.5group, drop=FALSE)+theme_light()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),axis.text = element_blank(),
          axis.ticks =  element_blank(),
          panel.background = element_rect(colour = "black", size=1))+ggtitle(TCR2plot$clonotype_id.consensus[i])
  pdf(paste0(TCR2plot$clonotype_id.consensus[i],".5groups.UMAP.pdf"),width=5,height=3)
  print(p1)
  dev.off()
  
}

