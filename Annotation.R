library(Seurat)
library(dplyr)
library(SingleR)
library(ggplot2)
library(tibble)
load("MFI_merge_singlet.Rdata")
#####SingleR#####
load("C4.Rdata")
C4_mtx <- GetAssayData(C4, slot="data")
C4_mtx <- log1p(C4_mtx)
C4$anno_label <- ifelse(C4$annotation%in%c("DC1","DC2"),"DC",
                        ifelse(C4$annotation%in%c("dM1","dM2","dM3"),"dM",
                            ifelse(C4$annotation%in%c("Granulocytes"),"Granulocytes",
                                ifelse(C4$annotation%in%c("HB"),"HB",
                                    ifelse(C4$annotation%in%c("ILC3"),"ILC3",
                                        ifelse(C4$annotation%in%c("MO"),"MO",
                                            ifelse(C4$annotation%in%c("Plasma"),"Plasma",
                                                ifelse(C4$annotation%in%c("Tcells"),"T cells",
                                                    ifelse(C4$annotation%in%c("dNK p","dNK1","dNK2","dNK3","NK CD16-","NK CD16+"),"NK",
                                                        ifelse(C4$annotation%in%c("dP1","dP2"),"dP",
                                                            ifelse(C4$annotation%in%c("dS1","dS2","dS3"),"dS",
                                                                ifelse(C4$annotation%in%c("Endo (f)","Endo (m)","Endo L"),"Endo",
                                                                    ifelse(C4$annotation%in%c("Epi1","Epi2"),"Epi",
                                                                        ifelse(C4$annotation%in%c("fFB1","fFB2"),"fFB",
                                                                            ifelse(C4$annotation%in%c("EVT"),"EVT",
                                                                                ifelse(C4$annotation%in%c("VCT"),"VCT",
                                                                                    ifelse(C4$annotation%in%c("SCT"),"SCT","")))))))))))))))))
singleR_res <- SingleR(test = log1p(GetAssayData(MFI_merge_singlet,assay = "RNA",slot = "data")),method = "cluster",clusters = Idents(MFI_merge_singlet), ref = C4_mtx, labels = C4$anno_label)
anno_df <- as.data.frame(singleR_res$labels)
anno_df$Index <- rownames(singleR_res)
colnames(anno_df)[1] <- 'SingleR_anno'
MFI_merge_singlet$SingleR_anno <- anno_df$SingleR_anno[match(Idents(MFI_merge_singlet), anno_df$Index)]

#####Annotation and subcluster-----
anno_df <- read.csv("./new_res/anno_df_28clus.csv",row.names = 1)
current.cluster.ids <- anno_df$Index
new.cluster.ids <- anno_df$Annotation
MFI_merge_singlet@active.ident <- plyr::mapvalues(x = MFI_merge_singlet@active.ident,
                                                  from = current.cluster.ids,
                                                  to = new.cluster.ids)
for(clus in c("DC","dP","dS","fFB","Endo","Epi","NK","T cells","dM","MO")){
  dir.create(paste0("./new_res/clus_",clus))
  MFI_subset <- subset(MFI_merge_singlet,idents = clus)
  MFI_subset <- RunPCA(MFI_subset, verbose = FALSE,assay = "integrated")
  MFI_subset <- RunUMAP(MFI_subset,assay = "integrated",dims = 1:20)
  MFI_subset <- RunTSNE(MFI_subset,assay = "integrated",dims = 1:20)
  MFI_subset <- FindNeighbors(MFI_subset, dims = 1:20)
  MFI_subset <- FindClusters(MFI_subset, verbose = FALSE,resolution = 0.1)
  save(MFI_subset,file = paste0("./new_res/clus_",clus,"/MFI_",clus,".Rdata"))
}
for(clus in c("DC","dP","dS","fFB","Endo","Epi","NK","T cells","dM","MO")){
  load(paste0("./new_res/clus_",clus,"/MFI_",clus,".Rdata"))
  MFI.subset.markers <- FindAllMarkers(object = MFI_subset,
                                       only.pos = TRUE,
                                       min.pct = 0.25,
                                       logfc.threshold = 0.25)
}
for(clus in c("DC","dP","dS","fFB","Endo","Epi","NK","T cells","dM","MO")){
  load(paste0("./new_res/clus_",clus,"/MFI_",clus,".Rdata"))
  MFI.subset.meta <- MFI_subset@meta.data
}

#####subset T cells-----
MFI_singlet_metadata <- read.csv("./new_res/MFI_singlet_metadata_0319_update.csv",row.names = 1)
MFI_merge_singlet <- AddMetaData(MFI_merge_singlet,MFI_singlet_metadata)
MFI_merge_singlet <- SetIdent(MFI_merge_singlet,cells = colnames(MFI_merge_singlet),value = MFI_merge_singlet$final_anno)
#T cells
MFI_subset <- subset(MFI_merge_singlet,idents = "Tcells")
MFI_subset <- RunPCA(MFI_subset, verbose = FALSE,assay = "integrated")
MFI_subset <- RunUMAP(MFI_subset,assay = "integrated",dims = 1:20)
MFI_subset <- RunTSNE(MFI_subset,assay = "integrated",dims = 1:20)
MFI_subset <- FindNeighbors(MFI_subset, dims = 1:20)
MFI_T <- FindClusters(MFI_subset, verbose = FALSE,resolution = 0.3)
save(MFI_T,file = "./new_res/MFI_T.Rdata")

#####FindMarkers#####
MFI.anno.markers <- FindAllMarkers(object = MFI_merge_singlet,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25)