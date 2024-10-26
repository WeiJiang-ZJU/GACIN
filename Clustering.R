library(Seurat)
library(dplyr)
library(ggplot2)
load("MFI_merge.Rdata")
#####Clustering#####
MFI_merge <- RunTSNE(MFI_merge, dims = 1:30)
MFI_merge <- FindNeighbors(MFI_merge, reduction = "pca", dims = 1:30)
MFI_merge <- FindClusters(MFI_merge, resolution = 0.6)

#####DoubletFinder#####
library(DoubletFinder)
load("./C0.Rdata")
C0[["percent.mt"]] <- PercentageFeatureSet(C0, pattern = "^MT-")
C0 <- subset(C0, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C1.Rdata")
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")
C1 <- subset(C1, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C2.Rdata")
C2[["percent.mt"]] <- PercentageFeatureSet(C2, pattern = "^MT-")
C2 <- subset(C2, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C3.Rdata")
C3[["percent.mt"]] <- PercentageFeatureSet(C3, pattern = "^MT-")
C3 <- subset(C3, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C4.Rdata")
load("./C5.Rdata")
C5[["percent.mt"]] <- PercentageFeatureSet(C5, pattern = "^MT-")
C5 <- subset(C5, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C6.Rdata")
C6[["percent.mt"]] <- PercentageFeatureSet(C6, pattern = "^MT-")
C6 <- subset(C6, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C7.Rdata")
C7[["percent.mt"]] <- PercentageFeatureSet(C7, pattern = "^MT-")
C7 <- subset(C7, subset = nFeature_RNA > 200 & percent.mt < 30)
load("./C8.Rdata")
C8[["percent.mt"]] <- PercentageFeatureSet(C8, pattern = "^MT-")
C8 <- subset(C8, subset = nFeature_RNA > 200 & percent.mt < 30)
{
    C8 <- NormalizeData(C8)
    C8 <- FindVariableFeatures(C8, selection.method = "vst", nfeatures = 2000)
    C8 <- ScaleData(C8, features = rownames(C8))
    C8 <- RunPCA(C8, npcs = 30, verbose = FALSE)
    C8 <- RunUMAP(C8, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C8, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C8)/7*8*1e-6
    C8$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C8),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C8$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(C8))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C8 <- doubletFinder_v3(C8, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C8_metadata <- C8@meta.data
    save(C8_metadata,file = "C8_metadata.Rdata")
}
{
    C7 <- NormalizeData(C7)
    C7 <- FindVariableFeatures(C7, selection.method = "vst", nfeatures = 2000)
    C7 <- ScaleData(C7, features = rownames(C7))
    C7 <- RunPCA(C7, npcs = 30, verbose = FALSE)
    C7 <- RunUMAP(C7, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C7, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C7)/4*8*1e-6
    C7$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C7),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C7$seurat_clusters)   
    nExp_poi <- round(DoubletRate*ncol(C7))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C7 <- doubletFinder_v3(C7, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C7_metadata <- C7@meta.data
    save(C7_metadata,file = "C7_metadata.Rdata")
}
{
    C6 <- NormalizeData(C6)
    C6 <- FindVariableFeatures(C6, selection.method = "vst", nfeatures = 2000)
    C6 <- ScaleData(C6, features = rownames(C6))
    C6 <- RunPCA(C6, npcs = 30, verbose = FALSE)
    C6 <- RunUMAP(C6, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C6, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C6)/3*8*1e-6
    C6$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C6),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C6$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(C6))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C6 <- doubletFinder_v3(C6, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C6_metadata <- C6@meta.data
    save(C6_metadata,file = "C6_metadata.Rdata")
}
{
    C5 <- NormalizeData(C5)
    C5 <- FindVariableFeatures(C5, selection.method = "vst", nfeatures = 2000)
    C5 <- ScaleData(C5, features = rownames(C5))
    C5 <- RunPCA(C5, npcs = 30, verbose = FALSE)
    C5 <- RunUMAP(C5, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C5, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C5)/3*8*1e-6
    C5$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C5),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C5$seurat_clusters)  
    nExp_poi <- round(DoubletRate*ncol(C5))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C5 <- doubletFinder_v3(C5, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C5_metadata <- C5@meta.data
    save(C5_metadata,file = "C5_metadata.Rdata")
}
{
    C4 <- NormalizeData(C4)
    C4 <- FindVariableFeatures(C4, selection.method = "vst", nfeatures = 2000)
    C4 <- ScaleData(C4, features = rownames(C4))
    C4 <- RunPCA(C4, npcs = 30, verbose = FALSE)
    C4 <- RunUMAP(C4, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C4, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C4)/7*8*1e-6
    C4$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C4),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C4$seurat_clusters)  
    nExp_poi <- round(DoubletRate*ncol(C4))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C4 <- doubletFinder_v3(C4, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C4_metadata <- C4@meta.data
    save(C4_metadata,file = "C4_metadata.Rdata")
}
{
    C3 <- NormalizeData(C3)
    C3 <- FindVariableFeatures(C3, selection.method = "vst", nfeatures = 2000)
    C3 <- ScaleData(C3, features = rownames(C3))
    C3 <- RunPCA(C3, npcs = 30, verbose = FALSE)
    C3 <- RunUMAP(C3, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C3, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C3)*8*1e-6
    C3$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C3),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C3$seurat_clusters)  
    nExp_poi <- round(DoubletRate*ncol(C3))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C3 <- doubletFinder_v3(C3, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C3_metadata <- C3@meta.data
    save(C3_metadata,file = "C3_metadata.Rdata")
}
{
    C2 <- NormalizeData(C2)
    C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures = 2000)
    C2 <- ScaleData(C2, features = rownames(C2))
    C2 <- RunPCA(C2, npcs = 30, verbose = FALSE)
    C2 <- RunUMAP(C2, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C2, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C2)/10*8*1e-6
    C2$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C2),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C2$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(C2))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C2 <- doubletFinder_v3(C2, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C2_metadata <- C2@meta.data
    save(C2_metadata,file = "C2_metadata.Rdata")
}
{
    C1 <- NormalizeData(C1)
    C1 <- FindVariableFeatures(C1, selection.method = "vst", nfeatures = 2000)
    C1 <- ScaleData(C1, features = rownames(C1))
    C1 <- RunPCA(C1, npcs = 30, verbose = FALSE)
    C1 <- RunUMAP(C1, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(C1, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(C1)/6*8*1e-6
    C1$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(C1),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(C1$seurat_clusters)   
    nExp_poi <- round(DoubletRate*ncol(C1))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    C1 <- doubletFinder_v3(C1, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    C1_metadata <- C1@meta.data
    save(C1_metadata,file = "C1_metadata.Rdata")
}
{
    C0_split <- SplitObject(C0,split.by = "donor")
    for (ID in names(C0_split)) {
    tmp <- C0_split[[ID]]
    tmp <- NormalizeData(tmp)
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    tmp <- ScaleData(tmp, features = rownames(tmp))
    tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
    tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:30)
    sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    DoubletRate = ncol(tmp)*8*1e-6
    tmp$seurat_clusters <- MFI_merge$seurat_clusters[match(colnames(tmp),colnames(MFI_merge))]
    homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
    nExp_poi <- round(DoubletRate*ncol(tmp))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pK_bcmvn,
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    tmp_metadata <- tmp@meta.data
    save(tmp_metadata,file = paste0(ID,"_metadata.Rdata"))
    }
}

#####DoubletFinder-merge#####
MFI_merge$Doublet_res <- ""
MFI_merge_metadata <- MFI_merge@meta.data
for (i in c(paste0("C0-P",c(1:37)))) {
  load(paste0(i,"_metadata.Rdata"))
  tmp <- tmp_metadata
  rownames(tmp) <- paste0(rownames(tmp),"_1")
  if(all(rownames(tmp)==rownames(MFI_merge_metadata[which(MFI_merge_metadata$donor==i),]))==TRUE)
  {
    ID <- data.frame(row.names = rownames(tmp),Doublet_res = tmp[,ncol(tmp)])
    MFI_colnames_subset <- colnames(MFI_merge)[which(colnames(MFI_merge)%in%rownames(ID))]
    MFI_merge$Doublet_res[which(colnames(MFI_merge)%in%rownames(ID))] <- ID$Doublet_res[match(MFI_colnames_subset,rownames(ID))]
  }else{print(paste0(i,"_False"))}
}
load("./MFI_merge_metadata.Rdata")
MFI_merge <- AddMetaData(MFI_merge,metadata = MFI_merge_metadata)
MFI_merge <- SetIdent(MFI_merge,cells = colnames(MFI_merge),value = MFI_merge$Doublet_res)
MFI_merge_singlet <- subset(MFI_merge,idents = "Singlet")
MFI_merge_singlet <- SetIdent(MFI_merge_singlet,cells = colnames(MFI_merge_singlet),value = MFI_merge_singlet$seurat_clusters)

#####FindMarkers#####
MFI.markers <- FindAllMarkers(object = MFI_merge_singlet,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)