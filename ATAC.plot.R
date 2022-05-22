

library(Signac)
library(Seurat)
library(GenomicRanges)

library(S4Vectors)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(harmony)

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(reshape2)
library(gridExtra)
color.5group<-c("#bea0d2", "#e6194b" ,"#8eb166", "#66b3dd", "#4363d8")

setwd("PATH_TO_ATAC_DATA")
sample.names<-list.dirs(getwd(),full.names=F,recursive=F)

####################
# combine peaks
####################

atac.peak.gr.list<-list()
for(i in 1:length(sample.names)){
  atac.peaks<-read.table(
    file = paste0(sample.names[i],"/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  atac.peak.gr.list[[i]] <- makeGRangesFromDataFrame(atac.peaks)
}
names(atac.peak.gr.list)<-sample.names
combined.peaks <- reduce(x = c(atac.peak.gr.list[[1]],atac.peak.gr.list[[2]],atac.peak.gr.list[[3]],atac.peak.gr.list[[4]],
                               atac.peak.gr.list[[5]],atac.peak.gr.list[[6]],atac.peak.gr.list[[7]],atac.peak.gr.list[[8]] ))
# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[ peakwidths > 20]
combined.peaks


md.list<-c()
for(i in 1:length(sample.names)){
  md.list[[i]] <- read.table(
    file = paste0(sample.names[i],"/per_barcode_metrics.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  md.list[[i]]<-md.list[[i]][which(md.list[[i]]$is_cell==1),]
}



####################
# create fragment objects
####################
atac.frag.list<-list()
atac.counts.list<-list()

#for(i in 1:length(sample.names)){
for(i in 1:length(sample.names)){ 
  atac.frag.list[[i]]<-CreateFragmentObject(
    path =paste0(sample.names[i],"/atac_fragments.tsv.gz") ,
    cells = rownames(md.list[[i]]))
  cells<-Cells(atac.frag.list[[i]])
  names(cells) <- cells
  Cells(atac.frag.list[[i]]) <- cells
  
  atac.counts.list[[i]] <- FeatureMatrix(
    fragments = atac.frag.list[[i]],
    features = combined.peaks,
    cells = rownames(md.list[[i]])
  )
}

########################
# Create the objects
########################
# add hg38 annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

atac_assay.list<-list()
atac.list<-list()

for(i in 1:length(sample.names)){
  atac_assay.list[[i]]<-CreateChromatinAssay(atac.counts.list[[i]], fragments = atac.frag.list[[i]])
  atac.list[[i]] <- CreateSeuratObject(atac_assay.list[[i]], assay = "ATAC", meta.data = md.list[[i]])
  atac.list[[i]]$dataset<-sample.names[i]
  atac.list[[i]]$patient<-substr(sample.names[i],1,7)
  atac.list[[i]]$pct_reads_in_peaks <-atac.list[[i]]$atac_peak_region_fragments / atac.list[[i]]$nCount_ATAC * 100
  Annotation(atac.list[[i]]) <- annotations
} 

# merge all samples
atac.combined <- merge(
  x = atac.list[[1]],
  y = list(atac.list[[2]],atac.list[[3]],atac.list[[4]],atac.list[[5]],atac.list[[6]],atac.list[[7]],atac.list[[8]]),
  add.cell.ids = sample.names
)

# add activity
gene.activities <- GeneActivity(atac.combined)
atac.combined[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
atac.combined <- NormalizeData(
  object = atac.combined,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize')



#basic QC
atac.combined <- NucleosomeSignal(object = atac.combined)
# compute TSS enrichment score per cell
atac.combined <- TSSEnrichment(object = atac.combined, fast = FALSE)
atac.combined$high.tss <- ifelse(atac.combined$TSS.enrichment > 2, 'High', 'Low')

# remove low quality cells
DefaultAssay(atac.combined) <- 'ATAC'
atac.4tumors.QC <- subset(
  x = atac.4tumors,
  subset = atac_peak_region_fragments > 3000 &
    atac_peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)


####################
# harmony integration
###################
atac.combined <- RunTFIDF(atac.combined)
atac.combined <- FindTopFeatures(atac.combined, min.cutoff = 20)
atac.combined <- RunSVD(atac.combined)
atac.combined <- RunUMAP(atac.combined, dims = 2:40, reduction = 'lsi')


library(harmony)
atac.4tumors.harmony.integr <- RunHarmony(
  object = atac.combined,
  group.by.vars = 'patient',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE,plot_convergence = TRUE
)

atac.4tumors.harmony.integr <- RunUMAP(atac.4tumors.harmony.integr, dims = 1:40, reduction = 'harmony')
p1 <- DimPlot(atac.4tumors.harmony.integr, group.by = 'dataset', pt.size = 0.1, cols=color3) + ggplot2::ggtitle("Harmony integration")
pdf("harmony.integrated.Sample.pdf")
print(p1)
dev.off()


############################################################################################################
# analyze TF motif using signac, JASPAR2020, and motifmatchr
############################################################################################################

library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library("BSgenome.Hsapiens.UCSC.hg38")
library(motifmatchr)
set.seed(1234)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
genome <- keepBSgenomeSequences(BSgenome.Hsapiens.UCSC.hg38, sequences_to_keep)


keep.peaks <- as.logical(as.character(seqnames(granges(atac.4tumors.harmony.integr)))%in%sequences_to_keep)
atac.4tumors.harmony.integr <- atac.4tumors.harmony.integr[keep.peaks, ]

atac.4tumors.harmony.integr<- AddMotifs(
  object = atac.4tumors.harmony.integr,
  genome = genome,
  pfm =  pfm
)

############################################################################################################
# add T cell differentiation states identified from paired scRNA-seq assay
############################################################################################################



############################################################################################################
# plot acticity scores for T differentiation states marker genes
############################################################################################################
atac.4tumors.QC.harmony.integr.subset.group$path<-factor(atac.4tumors.QC.harmony.integr.subset.group$path,levels=c("Stemlike_Recirc","Stemlike_TRM",
                                                                                                                     "Effector_TRM","Proliferation_TRM",
                                                                                                                     "Exhaustion_TRM"))
atac.4tumors.QC.harmony.integr.subset.group$Group<-mapvalues(atac.4tumors.QC.harmony.integr.subset.group$path,
                                                             from=levels(atac.4tumors.QC.harmony.integr.subset.group$path),
                                                             to=c("Stemlike_recirc","Stemlike_TRM","Effector_TRM","Proliferative_TRM","Exhausted_TRM"))
celltype.genelist<-list(Stemness=c("SELL","LEF1","CD28","CD27","CCR7","IL7R","CXCR5","TCF7","BACH2","JUNB","EGR1","KLF2"),
                        Effector=c("GZMK","GZMH","GZMB","PRF1","GNLY","IFNG","FASLG","FGFBP2"),
                        Proliferation=c("TOP2A","MKI67","CDK1","STMN1","DNMT1","MCM7"),
                        Exhaustion=c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","CD74","CXCL13","LAYN","MYO7A","HLA-DRB1","HLA-DQA1",
                                     "TOX","TOX2","BATF","ETV1","ID2","ZNF683","RBPJ","TBX21","RUNX1","RUNX3","RUNX2","EZH2"))



df2<- data.frame(ID = rep(names(celltype.genelist[1:length(celltype.genelist)]), sapply(celltype.genelist[1:length(celltype.genelist)], length)),
                 Obs = unlist(celltype.genelist[1:length(celltype.genelist)]))
#df<-rbind(df,df1)

gene_list <- rownames(atac.4tumors.QC.harmony.integr.subset.group@assays$ACTIVITY@data)
df2<- df2[df2$Obs %in% gene_list, ]
df2<- df2[!duplicated(df2$Obs),]
select_gene<- df2$Obs
gene.categories<-factor(df2$ID, levels=unique(df2$ID)) # define categories for genes
color.a<-c(color.5group)[as.numeric(gene.categories)] 

pos<-c()
d<-c()
for(i in 1:length(levels(gene.categories))){
  if(i==1){
    d=0
  }else{d= d+length(which(gene.categories==levels(gene.categories)[i-1]))}
  pos<-c(pos,d+length(which(gene.categories==levels(gene.categories)[i]))/2)
}



DefaultAssay(atac.4tumors.QC.harmony.integr.subset.group)<-"ACTIVITY"
gene.cat.lable<-data.frame(pos=pos, y=rep(length(unique(atac.4tumors.QC.harmony.integr.subset.group$Group))+1,length(pos)), cat=levels(gene.categories))
gene.cat.lable$cat<-factor(gene.cat.lable$cat)

p1<-DotPlot(atac.4tumors.QC.harmony.integr.subset.group, group.by = "Group", features=select_gene, )+ RotatedAxis()+scale_colour_gradient2(low ="#4575B4", mid = "#FFFFBF", high = "#D73027")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9, colour =color.a))+scale_y_discrete(expand = expansion(mult = c(0.05, .2)))+
  geom_text(data = gene.cat.lable, aes(x = pos, y = y, label = cat), colour=color3[c(1:length(celltype.genelist))],angle=30, vjust=3)




