library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
mtx <- read.csv("./cyToF/tsne.csv",row.names = 1,check.names = F)
colnames(mtx)
mtx_data <- t(mtx[,c(1:41)])
mtx_meta <- mtx[,c(42:44)]
cytof <- CreateSeuratObject(counts = mtx_data,meta.data = mtx_meta,min.cells = 0,min.features = 0)
#QC
VlnPlot(cytof,c("nCount_RNA","nFeature_RNA"),pt.size = 0)
cytof <- subset(cytof, subset = nFeature_RNA > 3 & nCount_RNA > 10)
cytof <- ScaleData(cytof, features = rownames(cytof),do.scale = F)
cytof <- RunPCA(cytof, features = rownames(cytof))
ElbowPlot(cytof,ndims = 30)
cytof <- FindNeighbors(cytof, dims = 1:11)
cytof <- FindClusters(cytof, resolution = 1)
cytof <- RunUMAP(cytof,dims = 1:11)
cytof <- RunTSNE(cytof,dims = 1:11)
current.cluster.ids <- levels(Idents(cytof))
new.cluster.ids <- c("Fibroblast","Trophoblast","Myeloid cells","T cells",
                     "Myeloid cells","T cells","Fibroblast","Epithelia",
                     "NK cells","Myeloid cells","Myeloid cells","Myeloid cells",
                     "Trophoblast","Fibroblast","Myeloid cells","Epithelia",
                     "Myeloid cells","Fibroblast","Myeloid cells","Trophoblast",
                     "Myeloid cells","Trophoblast","Endothelia","T cells",
                     "Endothelia","Fibroblast","Epithelia","Fibroblast",
                     "Endothelia")
cytof@active.ident <- plyr::mapvalues(x = cytof@active.ident, from = current.cluster.ids, 
                                      to = new.cluster.ids)
cytof$group <- gsub("G1","Early",cytof$group)
cytof$group <- gsub("G2","Middle",cytof$group)
cytof$group <- gsub("G3","Late",cytof$group)
cytof$group <- factor(cytof$group,levels = c("Early","Middle","Late"))

Fibro <- subset(cytof,idents = "Fibroblast")
Fibro <- RunPCA(Fibro, features = rownames(Fibro))
ElbowPlot(Fibro,ndims = 30)
Fibro <- FindNeighbors(Fibro, dims = 1:10)
Fibro <- FindClusters(Fibro, resolution = 0.5)
Fibro <- RunTSNE(Fibro,dims = 1:10)
highCells=colnames(subset(x = Fibro, subset = CD90 > 1, slot = 'counts'))
highORlow=ifelse(colnames(Fibro) %in% highCells,'high','low')
table(highORlow)
Fibro@meta.data$highORlow_CD90=highORlow
Fibro$subclus = as.character(Idents(Fibro))
Fibro$subclus[which(Idents(Fibro)=="8" & Fibro$highORlow_CD90=="low")] <- "11"
Fibro <- SetIdent(Fibro,cells = colnames(Fibro),value = Fibro$subclus)
current.cluster.ids <- as.character(c(0:11))
new.cluster.ids <- c("dP","fFB2","dP","dP","dP","dP","fFB1","dP","dS1","Fibroblast_1","Fibroblast_2","dS2")
Fibro@active.ident <- plyr::mapvalues(x = Fibro@active.ident, from = current.cluster.ids, 
                                      to = new.cluster.ids)
Idents(Fibro) <- factor(Idents(Fibro),levels = c("dS1","dS2","dP","fFB1","fFB2","Fibroblast_1","Fibroblast_2"))

Fibro$group <- gsub("G1","Early",Fibro$group)
Fibro$group <- gsub("G2","Middle",Fibro$group)
Fibro$group <- gsub("G3","Late",Fibro$group)
Fibro$group <- factor(Fibro$group,levels = c("Early","Middle","Late"))
Fibro_anno <- Idents(Fibro)
save(Fibro_anno,file = "./cyToF/Fibro/Fibro_anno_dP.Rdata")

Tropho <- subset(cytof,idents = "Trophoblast")
Tropho <- RunPCA(Tropho, features = rownames(Tropho))
ElbowPlot(Tropho,ndims = 30)
Tropho <- FindNeighbors(Tropho, dims = 1:9)
Tropho <- FindClusters(Tropho, resolution = 0.1)
Tropho <- RunTSNE(Tropho,dims = 1:9)
current.cluster.ids <- levels(Idents(Tropho))
new.cluster.ids <- c("SCT","VCT","EVT")
Tropho@active.ident <- plyr::mapvalues(x = Tropho@active.ident, from = current.cluster.ids, 
                                       to = new.cluster.ids)
Tropho_anno <- Idents(Tropho)
save(Tropho_anno,file = "./cyToF/Tropho/Tropho_anno.Rdata")

Myeloid <- subset(cytof,idents = "Myeloid cells")
Myeloid <- RunPCA(Myeloid, features = rownames(Myeloid))
ElbowPlot(Myeloid,ndims = 30)
Myeloid <- FindNeighbors(Myeloid, dims = 1:12)
Myeloid <- FindClusters(Myeloid, resolution = 0.4)
Myeloid <- RunTSNE(Myeloid,dims = 1:12)
current.cluster.ids <- levels(Idents(Myeloid))
new.cluster.ids <- c("dM1","dM3","dM2","dM3","HB","DC2","HB","dM3","DC1","dM2","MO","HB","dM1")
Myeloid@active.ident <- plyr::mapvalues(x = Myeloid@active.ident, from = current.cluster.ids, 
                                       to = new.cluster.ids)
Myeloid_anno <- Idents(Myeloid)
save(Myeloid_anno,file = "./cyToF/Myeloid/Myeloid_anno.Rdata")

Tcells <- subset(cytof,idents = "T cells")
Tcells <- RunPCA(Tcells, features = rownames(Tcells))
ElbowPlot(Tcells,ndims = 30)
Tcells <- FindNeighbors(Tcells, dims = 1:12)
Tcells <- FindClusters(Tcells, resolution = 2)
Tcells <- RunTSNE(Tcells,dims = 1:12)
ident23 <- colnames(Tcells)[which(Idents(Tcells)=="23")]
Tcells <- FindClusters(Tcells, resolution = 0.8)
Tcells <- SetIdent(Tcells,cells = ident23,value = "14")
current.cluster.ids <- as.character(c(0:14))
new.cluster.ids <- c("NK CD16-","CD8+ T","CD8+ T","CD4+ T","Plasma","NK CD16+","NK CD16+","T cells 1","CD4+ T","Plasma","NK CD16+","CD4+ T","NK CD16-","Plasma","Treg")
Tcells@active.ident <- plyr::mapvalues(x = Tcells@active.ident, from = current.cluster.ids, 
                                      to = new.cluster.ids)
Tcells_anno <- Idents(Tcells)
save(Tcells_anno,file = "./cyToF/Tcells/Tcells_anno.Rdata")

Epi <- subset(cytof,idents = "Epithelia")
Epi <- RunPCA(Epi, features = rownames(Epi))
ElbowPlot(Epi,ndims = 30)
Epi <- FindNeighbors(Epi, dims = 1:7)
Epi <- FindClusters(Epi, resolution = 0.1)
Epi <- RunTSNE(Epi,dims = 1:7)
current.cluster.ids <- levels(Idents(Epi))
new.cluster.ids <- c("Epi2","Epi1","Epi2")
Epi@active.ident <- plyr::mapvalues(x = Epi@active.ident, from = current.cluster.ids, 
                                       to = new.cluster.ids)
Epi_anno <- Idents(Epi)
save(Epi_anno,file = "./cyToF/Epi/Epi_anno.Rdata")

Endo <- subset(cytof,idents = "Endothelia")
Endo <- RunPCA(Endo, features = rownames(Endo))
ElbowPlot(Endo,ndims = 30)
Endo <- FindNeighbors(Endo, dims = 1:9)
Endo <- FindClusters(Endo, resolution = 0.1)
Endo <- RunTSNE(Endo,dims = 1:9)
current.cluster.ids <- as.character(c(0:4))
new.cluster.ids <- c("Endo(f)","Endo(m)","Endo L","Endo L","Endo(m)")
Endo@active.ident <- plyr::mapvalues(x = Endo@active.ident, from = current.cluster.ids, 
                                       to = new.cluster.ids)
Endo_anno <- Idents(Endo)
save(Endo_anno,file = "./cyToF/Endo/Endo_anno.Rdata")

NK <- subset(cytof,idents = "NK cells")
NK <- RunPCA(NK, features = rownames(NK))
ElbowPlot(NK,ndims = 30)
NK <- FindNeighbors(NK, dims = 1:6)
NK <- FindClusters(NK, resolution = 0.5)
NK <- RunTSNE(NK,dims = 1:6)
current.cluster.ids <- as.character(c(0:8))
new.cluster.ids <- c("dNK2","dNK p","dNK3","dNK p","dNK3","dNK1","NK cells","T cells 2","T cells 3")
NK@active.ident <- plyr::mapvalues(x = NK@active.ident, from = current.cluster.ids, 
                                       to = new.cluster.ids)
NK_anno <- Idents(NK)
save(NK_anno,file = "./cyToF/NK/NK_anno.Rdata")

#####Merge#####
load("Fibro/Fibro_anno_dP.Rdata")
load("Endo/Endo_anno.Rdata")
load("Epi/Epi_anno.Rdata")
load("Myeloid/Myeloid_anno.Rdata")
load("NK/NK_anno.Rdata")
load("Tcells/Tcells_anno.Rdata")
load("Tropho/Tropho_anno.Rdata")
cytof <- SetIdent(cytof,cells = names(Endo_anno),value = Endo_anno)
cytof <- SetIdent(cytof,cells = names(Epi_anno),value = Epi_anno)
cytof <- SetIdent(cytof,cells = names(Fibro_anno),value = Fibro_anno)
cytof <- SetIdent(cytof,cells = names(Myeloid_anno),value = Myeloid_anno)
cytof <- SetIdent(cytof,cells = names(NK_anno),value = NK_anno)
cytof <- SetIdent(cytof,cells = names(Tcells_anno),value = Tcells_anno)
cytof <- SetIdent(cytof,cells = names(Tropho_anno),value = Tropho_anno)
cytof <- SetIdent(cytof,cells = colnames(cytof)[which(Idents(cytof)%in%c("dP1","dP2"))],value = "dP")
cytof$main_clus <- ifelse(Idents(cytof)%in%c("dP1","dP2"),"Perivascular cells",
                          ifelse(Idents(cytof)%in%c("dS1","dS2"),"Stroma",
                                 ifelse(Idents(cytof)%in%c("Epi1","Epi2"),"Epithelia",
                                        ifelse(Idents(cytof)%in%c("EVT","VCT","SCT"),"Trophoblast",
                                               ifelse(Idents(cytof)%in%c("fFB1","fFB2","Fibroblast_1","Fibroblast_2"),"Fibroblast",
                                                      ifelse(Idents(cytof)%in%c("Endo(m)","Endo(f)","Endo L"),"Endothelia","Immune cells"))))))
prop <- as.data.frame(table(Idents(cytof),cytof$sample))
cytof_meta <- cytof@meta.data
cytof_meta <- cytof_meta[,c(4,6)]
cytof_meta <- unique(cytof_meta)
prop$group <- cytof_meta$group[match(prop$Var2,cytof_meta$sample)]
ident_meta <- as.data.frame(table(Idents(cytof),cytof$main_clus))
ident_meta <- ident_meta[which(ident_meta$Freq!=0),]
prop$main_clus <- ident_meta$Var2[match(prop$Var1,ident_meta$Var1)]
prop <- prop %>% group_by(main_clus,Var2) %>% mutate(prop = Freq/sum(Freq))
library(gg.gap)
library(ggpubr)
prop_imm <- prop[which(prop$main_clus=="Immune cells"),]
prop_oth <- prop[which(prop$main_clus!="Immune cells"),]
p1 <- ggboxplot(prop_imm, x = "Var1", y = "prop",
                color = "black",fill = "group",width = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#00A087FF","#E64B35FF")) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.text = element_text(size = rel(0.8), margin = margin()),
        panel.spacing = unit(3, "pt"))+
  theme(axis.title = element_blank())+
  stat_compare_means(aes(group = group),method = "wilcox",hide.ns = T,label = "p.signif"
  )
p1 <- gg.gap(p1,ylim = c(0,0.9),segments = c(0.019,0.02))
p1
p2 <- ggboxplot(prop_oth, x = "Var1", y = "prop",
                color = "black",fill = "group",width = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#00A087FF","#E64B35FF")) +
  facet_grid(.~main_clus,scales = "free",space = "free")+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.text = element_text(size = rel(0.8), margin = margin()),
        panel.spacing = unit(3, "pt"))+
  theme(axis.title = element_blank())
p2 <- p2 +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     hide.ns = T,label = "p.signif"
  )
p2
ggarrange(p1, p2, common.legend = TRUE,nrow = 2,heights = c(2,1))
ggsave("./cyToF/cyToF_prop_all.pdf",width = 12,height = 10)

cytof_split <- SplitObject(cytof,split.by = "ident")
av <- list()
for (i in levels(Idents(cytof))) {
  tmp <- GetAssayData(cytof_split[[i]])
  av[[i]] <- as.data.frame(rowMeans(tmp))
  colnames(av[[i]]) <- i
}
av <- Reduce(cbind,av)
head(av)
library(pheatmap)
pairs.breaks <- seq(0, 5, by=0.1)
pheatmap(t(av),border_color = NA,scale = "none",cellwidth = 8,cellheight = 8,
         fontsize = 8,
         treeheight_row = 10,treeheight_col = 10,color = my_palette,
         filename = "./cyToF/heatmap_all_avg_noscale.pdf")
cytof_meta <- cytof@meta.data
cytof_meta$celltype <- Idents(cytof)
cytof_meta <- cytof_meta[,c(4,6,9)]
rownames(cytof_meta) <- colnames(cytof)
write.csv(cytof_meta,file = "./cyToF/cytof_meta.csv")