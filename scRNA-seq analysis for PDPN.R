# Install the R packages ####
library(limma)
library(Seurat)
library(dplyr)
library(plyr)
library(magrittr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(monocle)
library(viridis)
library(ComplexHeatmap)
library(CellChat)
library(ggpubr)
library(GOplot)
library(grid)
library(gridExtra)
library(forcats)
library(lattice)
library(fgsea)
# Import data ####
#DM1
rt=read.table("DM1.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_DM1=avereps(data)
DM1 <- CreateSeuratObject(counts = data_DM1,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
DM1[["percent.mt"]] <- PercentageFeatureSet(object = DM1, pattern = "^Mt-")
# VlnPlot(DM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#DM2
rt=read.table("DM2.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_DM2=avereps(data)
DM2 <- CreateSeuratObject(counts = data_DM2,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
DM2[["percent.mt"]] <- PercentageFeatureSet(object = DM2, pattern = "^Mt-")
# VlnPlot(DM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#PDPN1
rt=read.table("PDPN1.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_PDPN1=avereps(data)
PDPN1 <- CreateSeuratObject(counts = data_PDPN1,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
PDPN1[["percent.mt"]] <- PercentageFeatureSet(object = PDPN1, pattern = "^Mt-")
# VlnPlot(PDPN1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#PDPN2
rt=read.table("PDPN2.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_PDPN2=avereps(data)
PDPN2 <- CreateSeuratObject(counts = data_PDPN2,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
PDPN2[["percent.mt"]] <- PercentageFeatureSet(object = PDPN2, pattern = "^Mt-")
# VlnPlot(PDPN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#PDPN3
rt=read.table("PDPN3.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_PDPN3=avereps(data)
PDPN3 <- CreateSeuratObject(counts = data_PDPN3,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
PDPN3[["percent.mt"]] <- PercentageFeatureSet(object = PDPN3, pattern = "^Mt-")
# VlnPlot(PDPN3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#PDPN4
rt=read.table("PDPN4.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_PDPN4=avereps(data)
PDPN4 <- CreateSeuratObject(counts = data_PDPN4,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
PDPN4[["percent.mt"]] <- PercentageFeatureSet(object = PDPN4, pattern = "^Mt-")
# VlnPlot(PDPN4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#N1
rt=read.table("N1.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_N1=avereps(data)
N1 <- CreateSeuratObject(counts = data_N1,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
N1[["percent.mt"]] <- PercentageFeatureSet(object = N1, pattern = "^Mt-")
# VlnPlot(N1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#N2
rt=read.table("N2.expression_matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data_N2=avereps(data)
N2 <- CreateSeuratObject(counts = data_N2,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
N2[["percent.mt"]] <- PercentageFeatureSet(object = N2, pattern = "^Mt-")
# VlnPlot(N2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Data integration and bench information ####
# 
drg.big<- merge(N1, y = c(N2, DM1,DM2,PDPN1,PDPN2,PDPN3,PDPN4),project = "DRG")
# Read bentch file
batch <- read.table("batch.txt",sep="\t",check.names=F)
batch <- as.matrix(batch)
drg.big@meta.data$batch <- batch
# write.table(rownames(drg.big@meta.data),file="state.txt",sep="\t",row.names=F,quote=F)
# Group information in meta 
state <- read.table("state.txt",sep="\t",check.names=F)
state <- as.matrix(state)
drg.big@meta.data$state <- state
drg.list <- SplitObject(drg.big, split.by = "batch")
#
pdf(file="Genes per cell.pdf",width=18,height=6)
plot1 <- VlnPlot(object = drg.big, features = c("nFeature_RNA"), ncol = 1,pt.size=0,group.by = "state",cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
  xlab("") + ylab("Genes per cell") + ggtitle("")+theme(legend.position = "none")
plot2 <- VlnPlot(object = drg.big, features = c("nCount_RNA"), ncol = 1,pt.size=0,group.by = "state",cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
  xlab("") + ylab("UMIs per cell") + ggtitle("")+theme(legend.position = "none")
plot3 <- VlnPlot(object = drg.big, features = c("percent.mt"), ncol = 1,pt.size=0,group.by = "state",cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
  xlab("") + ylab("Mitochondrial genes ratio") + ggtitle("")
CombinePlots(plots = list(plot1, plot2, plot3),ncol = 3)
dev.off()
#
drg.list <- lapply(X = drg.list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#
features <- SelectIntegrationFeatures(object.list = drg.list)
anchors <- FindIntegrationAnchors(object.list = drg.list, anchor.features = features)
drg <- IntegrateData(anchorset = anchors)
DefaultAssay(drg) <- "integrated"
#
all.genes <- rownames(drg)
drg <- ScaleData(drg, features = all.genes)
drg <- RunPCA(drg, features = VariableFeatures(object = drg))
#
drg <- JackStraw(drg, num.replicate = 100)
drg <- ScoreJackStraw(drg, dims = 1:20)
# JackStrawPlot(drg, dims = 1:15)
# ElbowPlot(drg,ndims = 50)
drg <- FindNeighbors(drg, dims = 1:30)
drg <- FindClusters(drg, resolution = 1.2)
drg <- RunUMAP(drg, dims = 1:30)
drg <- RunTSNE(drg, dims = 1:30)
DimPlot(drg, reduction = "tsne",label = T)

# Clusters and Markers
drg.markers <- FindAllMarkers(drg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(drg.markers,file="drg.markers.txt",sep="\t",row.names=T,quote=F)

# Annotion
drg0 <- drg
new.cluster.ids <- c("Neuron",
                     "SGC",
                     "Neuron", 
                     "Neuron", 
                     "Neuron", 
                     "Neuron",
                     "SGC",
                     "Neuron", 
                     "SGC",
                     "SGC",
                     "SGC",
                     "Neuron", 
                     "VEC",
                     "Neuron", 
                     "Schwann Cell",
                     "Neuron", 
                     "Neuron", 
                     "Microglia",
                     "Neuron", 
                     "VSMC",
                     "PSGC",
                     "Fibroblast")
names(new.cluster.ids) <- levels(drg0)
drg0 <- RenameIdents(drg0, new.cluster.ids)
pdf(file="cluster.pdf",width=8,height=8)
DimPlot(drg0, reduction = "tsne", label = TRUE, pt.size = 1,label.size = 8) + NoLegend()+
  theme(legend.position = "none", 
        axis.title = element_text(size = 25), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

# DRG Cell types and markers ####
pdf(file="Marker.pdf",width=8,height=15)
chart1 <- FeaturePlot(drg, reduction = "tsne",features = c("Fabp7"),pt.size = 1,cols = c("white","blue")) + NoLegend() +theme(
        axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
chart2 <- FeaturePlot(drg, reduction = "tsne",features = c("Top2a"),pt.size = 1,cols = c("white","blue")) + NoLegend()  +theme(legend.position = "none", 
                                                                                                                              axis.title = element_blank(), 
                                                                                                                              axis.text  = element_blank(),
                                                                                                                              axis.ticks = element_blank(),
                                                                                                                              axis.line = element_blank())
chart3 <- FeaturePlot(drg, reduction = "tsne",features = c("Mpz"),pt.size = 1,cols = c("white","blue")) + NoLegend()  +theme(legend.position = "none", 
                                                                                                                                 axis.title = element_blank(), 
                                                                                                                                 axis.text  = element_blank(),
                                                                                                                                 axis.ticks = element_blank(),
                                                                                                                                 axis.line = element_blank())
chart4 <- FeaturePlot(drg, reduction = "tsne",features = c("Cldn5"),pt.size = 1,cols = c("white","blue")) + NoLegend()  +theme(legend.position = "none", 
                                                                                                                                   axis.title = element_blank(), 
                                                                                                                                   axis.text  = element_blank(),
                                                                                                                                   axis.ticks = element_blank(),
                                                                                                                                   axis.line = element_blank())
chart5 <- FeaturePlot(drg, reduction = "tsne",features = c("Lyz2"),pt.size = 1,cols = c("white","blue")) + NoLegend()  +theme(legend.position = "none", 
                                                                                                                                  axis.title = element_blank(), 
                                                                                                                                  axis.text  = element_blank(),
                                                                                                                                  axis.ticks = element_blank(),
                                                                                                                                  axis.line = element_blank())
chart6 <- FeaturePlot(drg, reduction = "tsne",features = c("Tpm2"),pt.size = 1,cols = c("white","blue")) + NoLegend()  +theme(legend.position = "none", 
                                                                                                                                  axis.title = element_blank(), 
                                                                                                                                  axis.text  = element_blank(),
                                                                                                                                  axis.ticks = element_blank(),
                                                                                                                                  axis.line = element_blank())
chart7 <- FeaturePlot(drg, reduction = "tsne",features = c("Dcn"),pt.size = 1,cols = c("white","blue")) + NoLegend() +theme(legend.position = "none", 
                                                                                                                                 axis.title = element_blank(), 
                                                                                                                                 axis.text  = element_blank(),
                                                                                                                                 axis.ticks = element_blank(),
                                                                                                                                 axis.line = element_blank())
grid.newpage()  ###新建图表版面
pushViewport(viewport(layout = grid.layout(4,2))) ####将版面分成2*2矩阵
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(chart1, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图chart3
print(chart2, vp = vplayout(1,2))     ###将(2,1)的位置画图chart2          
print(chart3, vp = vplayout(2,1))    ###将（2,2)的位置画图chart1
print(chart4, vp = vplayout(2,2))    ###将（2,2)的位置画图chart1
print(chart5, vp = vplayout(3,1))    ###将（2,2)的位置画图chart1
print(chart6, vp = vplayout(3,2))    ###将（2,2)的位置画图chart1
print(chart7, vp = vplayout(4,1))    ###将（2,2)的位置画图chart1
dev.off()

pdf(file="Legend.pdf",width=10,height=8)
FeaturePlot(drg, reduction = "tsne",features = c("Fabp7"),pt.size = 0.25,cols = c("white","blue")) 
dev.off()


# Subset neurons ####
Neuron <- subset(drg,idents = c("0","2","3","4","5","7","11","13","15","16","18"))
Neuron <- RunTSNE(Neuron, dims = 1:30)
Neuron0 <- Neuron
new.cluster.ids <- c("4",
                     "7",
                     "2", 
                     "1", 
                     "5", 
                     "6",
                     "9",
                     "3", 
                     "10",
                     "8",
                     "11")
names(new.cluster.ids) <- levels(Neuron0)
Neuron0 <- RenameIdents(Neuron, new.cluster.ids)
Neuron0@active.ident = factor(Neuron0@active.ident, levels = c("11","10", "9","8","7","6", "5","4","3","2","1"))

pdf(file="neuron_cluster.pdf",width=9,height=8)
DimPlot(Neuron0, reduction = "tsne", label = TRUE, pt.size = 2,label.size = 12)+
  theme(legend.position = "none",
        axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()



# Markers in dot plots图####
features = c("Calca","Ntrk1","Tac1","P2rx3","Mrgprb4","Nefh","Sst","Nppb","Exoc1l","Piezo2","Trpv1","Asic3","Fxyd7","Atp1b1","S100b")
pdf(file="neuron_dot.pdf",width=6,height=6)
DotPlot(Neuron0, features = features,cols = c("lightgrey", "red"))+
  theme(legend.position = "top",
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

# Neuronal types ####
new.cluster.ids <- c("NP",
                     "PEP",
                     "NP", 
                     "NP", 
                     "NP", 
                     "NP",
                     "C-LTMR",
                     "SOM", 
                     "MAAC",
                     "TRPM8",
                     "Unidentified")
names(new.cluster.ids) <- levels(Neuron)
Neuron <- RenameIdents(Neuron, new.cluster.ids)
Neuron@active.ident = factor(Neuron@active.ident, levels = c("NP",
                                                             "PEP",
                                                             "SOM",
                                                             "C-LTMR",
                                                             "TRPM8",
                                                             "MAAC",
                                                             "Unidentified"))
# Reorder
# Neuron@active.ident = factor(Neuron@active.ident, levels = c("Unidentified",
#                                                              "MAAC",
#                                                              "TRPM8",
#                                                              "C-LTMR",
#                                                              "PEP",
#                                                              "NP6",
#                                                              "NP5",
#                                                              "NP4",
#                                                              "NP3",
#                                                              "NP2",
#                                                              "NP1"))

pdf(file="neuron_ann_tsne.pdf",width=9,height=8)
DimPlot(Neuron, reduction = "tsne", label = TRUE, pt.size = 2,label.size = 13)+
  theme(axis.title = element_text(size = 25), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

pdf(file="neuron_ann_umap.pdf",width=9,height=8)
DimPlot(Neuron, reduction = "umap", label = TRUE, pt.size = 2,label.size = 13)+
  theme(axis.title = element_text(size = 25), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

# neuron.markers <- FindAllMarkers(Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.table(neuron.markers,file="neuron.markers.txt",sep="\t",row.names=T,quote=F)

# MAAC <- subset(Neuron,idents = c("MAAC"))
# Idents(MAAC) <- "state"
# neuron.markers <- FindMarkers(MAAC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0,ident.1 = "Diabetic rats with MA", ident.2 = "Diabetic rats without MA")
# write.table(neuron.markers,file="MAAC_DEGs.txt",sep="\t",row.names=T,quote=F)
# neuron.markers <- FindMarkers(MAAC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0,ident.1 = "Diabetic rats with MA", ident.2 = "Control group")
# write.table(neuron.markers,file="MAAC_DEGs2.txt",sep="\t",row.names=T,quote=F)

# PEP <- subset(Neuron,idents = c("PEP"))
# Idents(PEP) <- "state"
# neuron.markers <- FindMarkers(PEP, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0,ident.1 = "Diabetic rats with MA", ident.2 = "Diabetic rats without MA")
# write.table(neuron.markers,file="PEP_DEGs_PDPN_vs_DM.txt",sep="\t",row.names=T,quote=F)
# neuron.markers <- FindMarkers(PEP, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0,ident.1 = "Diabetic rats with MA", ident.2 = "Control group")
# write.table(neuron.markers,file="PEP_DEGs_PDPN_vs_N.txt",sep="\t",row.names=T,quote=F)

# Neuron markers in dot plot ####
features = c("Calca","Ntrk1","Tac1","P2rx3","Mrgprb4","Nefh","Sst","Nppb","Exoc1l","Piezo2","Trpv1","Asic3")
pdf(file="neuron_dot.pdf",width=6,height=6)
DotPlot(Neuron, features = features,cols = c("lightgrey", "red"))+
  theme(legend.position = "top",
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

# tSNE plot ####
# Change idents
Neuron@meta.data$cluster <- Neuron@active.ident
Idents(Neuron) <- Neuron$state
CG <- subset(Neuron,idents = c("Control group"))
pdf(file="control_group.pdf",width=9,height=9)
DimPlot(CG, reduction = "tsne", label = TRUE, pt.size = 3,label.size = 12,group.by = "cluster")+NoLegend()+
  theme(title = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()
pdf(file="control_group_gray.pdf",width=9,height=9)
DimPlot(CG, reduction = "tsne", label = F, pt.size = 3,label.size = 8,cols = "gray",group.by = "state")+NoLegend()+
  theme(axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()
WOMA <- subset(Neuron,idents = c("Diabetic rats without MA"))
pdf(file="Diabetic rats without MA.pdf",width=9,height=9)
DimPlot(WOMA, reduction = "tsne", label = F, pt.size = 3,label.size = 8,group.by = "cluster")+NoLegend()+
  theme(title = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()
WMA <- subset(Neuron,idents = c("Diabetic rats with MA"))
pdf(file="Diabetic rats with MA.pdf",width=9,height=9)
DimPlot(WMA, reduction = "tsne", label = F, pt.size = 3,label.size = 8,group.by = "cluster")+NoLegend()+
  theme(title = element_blank(),
        axis.title = element_blank(), 
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
dev.off()

# DEG heatmap ####
logFCfilter=0.25
adjPvalFilter=0.05
sig.markers=neuron.markers[(as.numeric(as.vector(neuron.markers$avg_log2FC))>logFCfilter & as.numeric(as.vector(neuron.markers$p_val_adj))<adjPvalFilter),]
top5 <- sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
maacmarkers <-sig.markers[sig.markers$cluster == "MAAC",]
maacmarkers = maacmarkers[order(maacmarkers$avg_log2FC,decreasing = TRUE),]
pdf(file="Heatmap_maac.pdf",width=10,height=20)
DoHeatmap(object = Neuron,label=T,
          features = maacmarkers$gene,
          size = 3.6, angle = -50, hjust=1,
          draw.lines=T) +
  scale_fill_gradientn(colors = c("skyblue", "white", "red"))+
  theme(legend.position = "top")
dev.off()

# DAVID GO enrichment ####
GO=read.table("GOBP.txt",sep="\t",header=T,check.names=F)
pdf(file="GO.pdf",width=12,height=6)
GO %>%
  mutate(Term = fct_reorder(Term,-log10(FDR))) %>%
  ggplot(aes(x=Term, y=-log10(FDR), fill = Fold_Enrichment)) +
  geom_bar(stat="identity",width = 0.5) +
  scale_fill_viridis_c(rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x>0.3, 
           scales::rescale(x,
                           to = to,
                           from = c(max(x, na.rm = TRUE), 0.3)),
           1)}) +
  theme_bw()+
  xlab(" ") +
  ylab("-log10(FDR)") +
  theme(legend.title = element_text(colour="black", size = 20))+
  theme(legend.text = element_text(size = 15,face = "plain"))+
  theme(axis.text.x = element_text(size=12,colour="black",angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(size=20,colour="black"),
        axis.title=element_text(size=20))+
  theme(plot.margin = unit(c(1,1,1,5), "cm"))
dev.off()

# correlation heatmap####
#
neuron.markers <- FindAllMarkers(Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top50 <- neuron.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
cluster.averages <- AverageExpression(Neuron,return.seurat = TRUE)
averages_matrix <- as.matrix(cluster.averages@assays$RNA@data)
select_matrix <-  averages_matrix[top50$gene,]
#
library(corrplot)
library("wesanderson")
pdf("corHeatmap.pdf",height=5,width=5)              
corrplot(corr=cor(select_matrix),
         method = "color",
         order = "FPC",
         tl.col="black",
         addCoef.col = "white",
         number.cex = 0.8,
         col= wes_palette("Zissou1", 100, type = "continuous"),
         col.lim = c(0.5,1),
         is.corr = F,
         addgrid.col = 'white'
)
dev.off()

# Cell count######
# table(Idents(Neuron))
# table(Neuron$orig.ident)
# table(Idents(drg),drg$orig.ident)
# table(Neuron$state)
# prop.table(table(Idents(Neuron),Neuron$state))
# Neuron$state = factor(Neuron$state, levels = c("Control group","Diabetic rats without MA","Diabetic rats with MA"))

# SCENIC ####
# loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
# library(SCopeLoomR)
# loom <- open_loom(loomPath)
# exprMat <- as.matrix(Neuron@assays$RNA@data)  #exprMat <- get_dgem(loom)
# cellInfo <- Neuron@active.ident #cellInfo <- get_cell_annotation(loom)
#close_loom(loom)
# dim(exprMat)
# head(cellInfo) 
cellInfo <- data.frame(cellInfo)
colnames(cellInfo) <-c("CellType")
# cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
# colVars <- list(CellType=c("microglia"="forestgreen", 
#                            "endothelial-mural"="darkorange", 
#                            "astrocytes_ependymal"="magenta4", 
#                            "oligodendrocytes"="hotpink", 
#                            "interneurons"="red3", 
#                            "pyramidal CA1"="skyblue", 
#                            "pyramidal SS"="darkblue"))
colVars <- list(CellType=c("NP1"="red3",
                           "NP2"="brown",
                           "NP3"="darkorange",
                           "NP4"="forestgreen",
                           "NP5"="darkgreen",
                           "NP6"="green2",
                           "PEP"="skyblue",
                           "C-LTMR"="darkblue",
                           "TRPM8"="purple",
                           "MAAC"="pink",
                           "Unidentified"="hotpink"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

library(SCENIC)
org <- "mgi" # or hgnc, or dmel
dbDir <- "~/Desktop/END" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on DRG" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
rm(exprMat)

runCorrelation(exprMat_filtered, scenicOptions)

#
## If launched in a new session, you will need to reload...
# setwd("...")
# loomPath <- "..."
# loom <- open_loom(loomPath)
# exprMat <- get_dgem(loom)
# close_loom(loom)
# genesKept <- loadInt(scenicOptions, "genesKept")
# exprMat_filtered <- exprMat[genesKept,]
# library(SCENIC)
# scenicOptions <- readRDS("int/scenicOptions.Rds")

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
close_loom(loom)
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
exprMat_log <- log2(exprMat+1)
dim(exprMat)

library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
#
exprMat_log <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)  
#
print(tsneFileName(scenicOptions))
#
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Dlx5", "Sox10", "Sox9","Irf1", "Stat6")],], plots="Expression")
# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()
#
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c( "Dlx5","Sox10")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)

regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
#
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Dlx5", "Irf1")]
#
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
#
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 
#
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Dlx5"]
viewMotifs(tableSubset) 
#
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
#
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

# Monocle analysis ####
mono_data <- subset(Neuron, idents = c('PEP',"MAAC"))
mono_data@meta.data$cluster <- mono_data@active.ident
#
matrix=as.matrix(mono_data@assays$RNA@data)
matrix=cbind(id=row.names(matrix),matrix)
write.table(matrix,file="Matrix.txt",quote=F,sep="\t",row.names=F)
#
sample=as.matrix(mono_data@meta.data)
sample=cbind(id=row.names(sample),sample)
write.table(sample,file="Sample.txt",quote=F,sep="\t",row.names=F)
#
geneAnn=data.frame(gene_short_name = row.names(matrix), row.names = row.names(matrix))
geneAnn=cbind(id=row.names(geneAnn),geneAnn)
write.table(geneAnn,file="Gene.txt",quote=F,sep="\t",row.names=F)
#
library(monocle)
expr_matrix=read.table("Matrix.txt",sep="\t",header=T,row.names=1,check.names=F)
sample_sheet=read.table("Sample.txt",sep="\t",header=T,row.names=1,check.names=F)
gene_annotation=read.table("Gene.txt",sep="\t",header=T,row.names=1,check.names=F)
# marker=read.table("Markers.txt", sep="\t",header=T,check.names=F)
expr_matrix = as.matrix(expr_matrix)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# 
DEG_genes <- differentialGeneTest(cds,fullModelFormulaStr = '~cluster',cores = 1)
# View
ordering_genes <- row.names(DEG_genes)[order(DEG_genes$pval)][1:153]
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <-reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
# cds <- orderCells(cds, root_state = 3);   
#
pdf(file="trajectory.pdf",width=8.5,height=6.5)
plot_cell_trajectory(cds, 
                     cell_size = 5,
                     color_by = "cluster",
                     show_branch_points = F,
                     show_state_number = F,
                     show_tree = T)+ 
  scale_color_manual(name = " ",breaks = c("PEP","PDPNAC"), values = c("#619CFF", "#F564E3"))+
  theme_classic(base_size = 16)+
  theme(legend.title = element_text(colour="black", size = 20))+
  theme(legend.text = element_text(colour="black", size = 20))+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  theme(legend.position = "top")
dev.off()

plot_cell_trajectory(cds, 
                     cell_size = 5,
                     color_by = "cluster",
                     show_branch_points = T,
                     show_state_number = T,
                     show_tree = T)
#
pdf(file="state_trajectory.pdf",width=8.5,height=6.5)
plot_cell_trajectory(cds, 
                     cell_size = 5,
                     color_by = "state",
                     show_branch_points = F,
                     show_state_number = F,
                     show_tree = T)+ 
  scale_color_manual(name = " ",breaks =  c("Control group","Diabetic rats without MA","Diabetic rats with MA"), 
                     values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_classic(base_size = 16)+
  theme(legend.title = element_text(colour="black", size = 20))+
  theme(legend.text = element_text(colour="black", size = 20))+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  theme(legend.position = "top")
dev.off()
#
pdf(file="Pseudotime.pdf",width=9,height=6.5)
plot_cell_trajectory(cds,cell_size = 5, color_by = "Pseudotime",show_branch_points = F)+
  theme(legend.position = "right")+
  theme_classic(base_size = 16)+
  theme(legend.title = element_text(colour="black", size = 20))+
  theme(legend.text = element_text(colour="black", size = 20))+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  theme(legend.position = "top")
dev.off()
# Heatmap
DEG<- row.names(subset(DEG_genes, qval < 0.3))
pdf(file="mono_heat2.pdf",width=8,height=20) 
plot_pseudotime_heatmap(cds[DEG,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = F)
dev.off()
# Output genes
p=plot_pseudotime_heatmap(cds[DEG,],
                          num_clusters = 4,
                          cores = 1,return_heatmap=T,
                          show_rownames = T)

clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.table(clustering,file="monodeg.txt",sep="\t",row.names=T,quote=F)



# CellChat ####
library(CellChat)
# Need Cairo for complexheatmap, then CellChat
# devtools::install_local("~/Downloads/CellChat-master.zip")
drgall <- drg
new.cluster.ids <- c("NP",
                     "SGC",
                     "PEP", 
                     "NP", 
                     "NP", 
                     "NP",
                     "SGC",
                     "NP", 
                     "SGC",
                     "SGC",
                     "SGC",
                     "C-LTMR", 
                     "VEC",
                     "SOM", 
                     "SC",
                     "MAAC", 
                     "TRPM8", 
                     "Microglia",
                     "Unidentified", 
                     "VSMC",
                     "PSGC",
                     "Fibroblast")
names(new.cluster.ids) <- levels(drgall)
drgall <- RenameIdents(drgall, new.cluster.ids)
# pdf(file="clusterall.pdf",width=8,height=8)
# DimPlot(drgall, reduction = "tsne", label = TRUE, pt.size = 1,label.size = 8) + NoLegend()+
#   theme(legend.position = "none",
#         axis.title = element_text(size = 25),
#         axis.text  = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank())
# dev.off()

#meta
drgall@meta.data$group <- drgall@active.ident
meta = drgall@meta.data
data.input = drgall@assays$RNA@data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "group") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#
df.net <- subsetCommunication(cellchat)
slot.name = "netP"
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
#
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#
groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
pdf(file="circle.pdf",width=6,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#
pdf(file="circle2.pdf",width=6,height=6)
par(mfrow = c(1,1), xpd=TRUE)
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[8, ] <- mat[8, ]
mat2[, 8] <- mat[, 8]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[8])
dev.off()

cellchat@netP$pathways

pathways.show <- c("PTN","MK","CALCR","ANGPTL","NT","SPP1","GALECTIN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#> [[1]]
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
pdf(file="bubble_source.pdf",width=5,height=5)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(2,5,10,12,13), remove.isolate = FALSE)
dev.off()
pdf(file="bubble_target.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use = c(2,5,13), targets.use = 8, remove.isolate = FALSE)
dev.off()
# violin plot
plotGeneExpression(cellchat, signaling = "PTN")
plotGeneExpression(cellchat, signaling = "CALCR")
plotGeneExpression(cellchat, signaling = "NT")

# network
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#
# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN"))
# #> Signaling role analysis on the cell-cell communication network from user's input
# gg1 + gg2
#
pdf(file="scatter.pdf",width=7,height=7)
netAnalysis_signalingRole_scatter(cellchat,font.size = 20,font.size.title = 20,label.size = 5,dot.size = c(5,8))
dev.off()
#
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2




