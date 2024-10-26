library(Seurat)
library(dplyr)
library(ggplot2)
#####DataSet0#####
meta <- read.csv("./data/DataSet0/new_metadata.csv",row.names = 1)
C0 <- Read10X(paste0("./data/DataSet0/",meta$sample[1],"/"))
C0 <- CreateSeuratObject(counts = C0,
                         min.cells = 3,
                         min.features = 200,
                         project = "Cohort0")
C0[["donor"]] <- "C0-P1"
C0[["tissue"]] <- meta$tissue[1]
C0[["group"]] <- "control"
C0[["cohort"]] <- "Cohort0"
C0[["pregnancy"]] <- meta$time[1]
for (ID in 2:37) {
  sampleID <- meta$sample[ID]
  tmp <- Read10X(paste0("./data/DataSet0/",sampleID,"/"))
  tmp <- CreateSeuratObject(counts = tmp,
                            min.cells = 3,
                            min.features = 200,
                            project = "Cohort0")
  tmp[["donor"]] <- paste0("C0-P",ID)
  tmp[["tissue"]] <- meta$tissue[ID]
  tmp[["group"]] <- "control"
  tmp[["cohort"]] <- "Cohort0"
  tmp[["pregnancy"]] <- meta$time[ID]
  C0 <- merge(C0,tmp)
}
save(C0,file = "C0.Rdata")

#####DataSet1#####
PD1 <- read.table("./data/DataSet1/D-20210705_matrix.tsv")
PP1 <- read.table("./data/DataSet1/P-20210705_matrix.tsv")
PD2 <- read.table("./data/DataSet1/E-20210713_matrix.tsv")
PP2 <- read.table("./data/DataSet1/P-20210713_matrix.tsv")
PD3 <- read.table("./data/DataSet1/D-20210805_matrix.tsv")
PP3 <- read.table("./data/DataSet1/P-20210805_matrix.tsv")
PD1 <- CreateSeuratObject(counts = PD1,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PP1 <- CreateSeuratObject(counts = PP1,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PD2 <- CreateSeuratObject(counts = PD2,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PP2 <- CreateSeuratObject(counts = PP2,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PD3 <- CreateSeuratObject(counts = PD3,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PP3 <- CreateSeuratObject(counts = PP3,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort1")
PD1[["donor"]] <- "C1-P1"
PP1[["donor"]] <- "C1-P1"
PD2[["donor"]] <- "C1-P2"
PP2[["donor"]] <- "C1-P2"
PD3[["donor"]] <- "C1-P3"
PP3[["donor"]] <- "C1-P3"
PD1[["tissue"]] <- "decidua"
PP1[["tissue"]] <- "placenta"
PD2[["tissue"]] <- "decidua"
PP2[["tissue"]] <- "placenta"
PD3[["tissue"]] <- "decidua"
PP3[["tissue"]] <- "placenta"
PD1[["group"]] <- "preeclampsia"
PP1[["group"]] <- "preeclampsia"
PD2[["group"]] <- "preeclampsia"
PP2[["group"]] <- "preeclampsia"
PD3[["group"]] <- "preeclampsia"
PP3[["group"]] <- "preeclampsia"
PD1[["cohort"]] <- "Cohort1"
PP1[["cohort"]] <- "Cohort1"
PD2[["cohort"]] <- "Cohort1"
PP2[["cohort"]] <- "Cohort1"
PD3[["cohort"]] <- "Cohort1"
PP3[["cohort"]] <- "Cohort1"
PD1[["pregnancy"]] <- "Late"
PP1[["pregnancy"]] <- "Late"
PD2[["pregnancy"]] <- "Late"
PP2[["pregnancy"]] <- "Late"
PD3[["pregnancy"]] <- "Late"
PP3[["pregnancy"]] <- "Late"
C1 <- merge(PD1,c(PP1,PD2,PP2,PD3,PP3))
C1
save(C1,file = "C1.Rdata")

#####DataSet2#####
FC1 <- Read10X("./data/DataSet2/C1/filtered_feature_bc_matrix/")
FC2 <- Read10X("./data/DataSet2/C2/filtered_feature_bc_matrix/")
FC3 <- Read10X("./data/DataSet2/C3/filtered_feature_bc_matrix/")
FC4 <- Read10X("./data/DataSet2/C4/filtered_feature_bc_matrix/")
FC5 <- Read10X("./data/DataSet2/C5/filtered_feature_bc_matrix/")
FP1 <- Read10X("./data/DataSet2/P1/filtered_feature_bc_matrix/")
FP2 <- Read10X("./data/DataSet2/P2/filtered_feature_bc_matrix/")
FP3 <- Read10X("./data/DataSet2/P3/filtered_feature_bc_matrix/")
FP4 <- Read10X("./data/DataSet2/P4/filtered_feature_bc_matrix/")
FP5 <- Read10X("./data/DataSet2/P5/filtered_feature_bc_matrix/")
FC1 <- CreateSeuratObject(counts = FC1,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FC2 <- CreateSeuratObject(counts = FC2,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FC3 <- CreateSeuratObject(counts = FC3,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FC4 <- CreateSeuratObject(counts = FC4,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FC5 <- CreateSeuratObject(counts = FC5,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FP1 <- CreateSeuratObject(counts = FP1,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FP2 <- CreateSeuratObject(counts = FP2,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FP3 <- CreateSeuratObject(counts = FP3,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FP4 <- CreateSeuratObject(counts = FP4,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FP5 <- CreateSeuratObject(counts = FP5,
                          min.cells = 3,
                          min.features = 200,
                          project = "Cohort2")
FC1[["donor"]] <- "C2-P1"
FC2[["donor"]] <- "C2-P2"
FC3[["donor"]] <- "C2-P3"
FC4[["donor"]] <- "C2-P4"
FC5[["donor"]] <- "C2-P5"
FP1[["donor"]] <- "C2-P6"
FP2[["donor"]] <- "C2-P7"
FP3[["donor"]] <- "C2-P8"
FP4[["donor"]] <- "C2-P9"
FP5[["donor"]] <- "C2-P10"
FC1[["tissue"]] <- "placenta"
FC2[["tissue"]] <- "placenta"
FC3[["tissue"]] <- "placenta"
FC4[["tissue"]] <- "placenta"
FC5[["tissue"]] <- "placenta"
FP1[["tissue"]] <- "placenta"
FP2[["tissue"]] <- "placenta"
FP3[["tissue"]] <- "placenta"
FP4[["tissue"]] <- "placenta"
FP5[["tissue"]] <- "placenta"
FC1[["group"]] <- "control"
FC2[["group"]] <- "control"
FC3[["group"]] <- "control"
FC4[["group"]] <- "control"
FC5[["group"]] <- "control"
FP1[["group"]] <- "preeclampsia"
FP2[["group"]] <- "preeclampsia"
FP3[["group"]] <- "preeclampsia"
FP4[["group"]] <- "preeclampsia"
FP5[["group"]] <- "preeclampsia"
FC1[["cohort"]] <- "Cohort2"
FC2[["cohort"]] <- "Cohort2"
FC3[["cohort"]] <- "Cohort2"
FC4[["cohort"]] <- "Cohort2"
FC5[["cohort"]] <- "Cohort2"
FP1[["cohort"]] <- "Cohort2"
FP2[["cohort"]] <- "Cohort2"
FP3[["cohort"]] <- "Cohort2"
FP4[["cohort"]] <- "Cohort2"
FP5[["cohort"]] <- "Cohort2"
C2 <- merge(FC1,c(FC2,FC3,FC4,FC5,FP1,FP2,FP3,FP4,FP5))
C2[["pregnancy"]] <- "Late"
C2
save(C2,file = "C2.Rdata")

#####DataSet3#####
MD <- read.table("./data/DataSet3/疾病组-GD01.expression_matrix.txt")
MC <- read.table("./data/DataSet3/正常组-ZD01.expression_matrix.txt")
MD <- CreateSeuratObject(counts = MD,
                         min.cells = 3,
                         min.features = 200,
                         project = "Cohort3")
MC <- CreateSeuratObject(counts = MC,
                         min.cells = 3,
                         min.features = 200,
                         project = "Cohort3")
MD[["donor"]] <- NA
MC[["donor"]] <- NA
MD[["tissue"]] <- "placenta"
MC[["tissue"]] <- "placenta"
MD[["group"]] <- "preeclampsia"
MC[["group"]] <- "control"
MD[["cohort"]] <- "Cohort3"
MC[["cohort"]] <- "Cohort3"
MD[["pregnancy"]] <- "Late"
MC[["pregnancy"]] <- "Late"
C3 <- merge(MD,MC)
C3
save(C3,file = "C3.Rdata")

#####DataSet4#####
C4 <- read.table("./data/DataSet4/raw_data_10x.txt",row.names = 1,header = T)
C4 <- C4[order(rowSums(C4),decreasing = T),]
genename_p <- unlist(lapply(strsplit(rownames(C4),split = "_"), function(X) X[[1]]))
C4 <- C4[!duplicated(genename_p),]
rownames(C4) <- unlist(lapply(strsplit(rownames(C4),split = "_"), function(X) X[[1]]))
meta <- read.table("./data/DataSet4/meta_10x.txt",sep = '\t')
meta <- meta[,c(1,2,4)]
colnames(meta) <- c("donor","tissue","annotation")
meta$tissue <- gsub("Decidua","decidua",meta$tissue)
meta$tissue <- gsub("Placenta","placenta",meta$tissue)
meta$tissue <- gsub("Blood","blood",meta$tissue)
C4 <- CreateSeuratObject(counts = C4,
                         min.cells = 3,
                         min.features = 200,
                         meta.data = meta,
                         project = "Cohort4")
C4[["group"]] <- "control"
C4[["cohort"]] <- "Cohort4"
C4[["pregnancy"]] <- "Early"
save(C4,file = "C4.Rdata")

#####DataSet5#####
load("./data/GSE130560/GSE130560_matrix.RData")
meta5 <- read.csv("./data/GSE130560/GSE130560_phenotype.csv",row.names = 1)
meta5 <- meta5[which(meta5$disease=="CTRL"),]
meta5 <- meta5[,c(5,9)]
colnames(meta5) <- c("donor","annotation_C5")
C5 <- matrix[,which(colnames(matrix)%in%rownames(meta5))]
C5 <- CreateSeuratObject(counts = C5,
                         min.cells = 3,
                         min.features = 200,
                         meta.data = meta5,
                         project = "Cohort5")
C5[["group"]] <- "control"
C5[["cohort"]] <- "Cohort5"
C5[["tissue"]] <- "decidua"
C5[["pregnancy"]] <- "Early"
save(C5,file = "C5.Rdata")

#####DataSet6#####
C6 <- read.table("./data/GSE171381/GSE171381_counts.tsv",sep = '\t',row.names = 1,header = T)
meta6 <- read.table("./data/GSE171381/GSE171381_metadata.tsv",row.names = 1,header = T)
meta6 <- meta6[which(meta6$covid == "cntrl"),]
meta6 <- meta6[,c(11,12,14)]
colnames(meta6) <- c("tissue","sample","annotation_C6")
C6 <- C6[,which(colnames(C6)%in%rownames(meta6))]
C6 <- CreateSeuratObject(counts = C6,
                         min.cells = 3,
                         min.features = 200,
                         meta.data = meta6,
                         project = "Cohort6")
C6[["group"]] <- "control"
C6[["cohort"]] <- "Cohort6"
C6[["pregnancy"]] <- "Late"
save(C6,file = "C6.Rdata")

#####DataSet7-----
GSM846 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945846/")
GSM847 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945847/")
GSM848 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945848/")
GSM849 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945849/")
GSM850 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945850/")
GSM851 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945851/")
GSM852 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945852/")
GSM853 <- Read10X("./data/DataSet7/GSE198373_RAW/GSM5945853/")
Middle1 <- CreateSeuratObject(counts = GSM846,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle2 <- CreateSeuratObject(counts = GSM847,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle3 <- CreateSeuratObject(counts = GSM848,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle4 <- CreateSeuratObject(counts = GSM849,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle5 <- CreateSeuratObject(counts = GSM850,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle6 <- CreateSeuratObject(counts = GSM851,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle7 <- CreateSeuratObject(counts = GSM852,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
Middle8 <- CreateSeuratObject(counts = GSM853,
                              min.cells = 3,
                              min.features = 200,
                              project = "Cohort7")
{
  Middle1[["donor"]] <- "C7-P1"
  Middle1[["tissue"]] <- "placenta"
  Middle1[["cohort"]] <- "Cohort7"
  Middle1[["group"]] <- "control"
  Middle1[["pregnancy_week"]] <- 24
  Middle1[["position"]] <- "VC"
  Middle2[["donor"]] <- "C7-P2"
  Middle2[["tissue"]] <- "placenta"
  Middle2[["cohort"]] <- "Cohort7"
  Middle2[["group"]] <- "control"
  Middle2[["pregnancy_week"]] <- 24
  Middle2[["position"]] <- "SC"
  Middle3[["donor"]] <- "C7-P3"
  Middle3[["tissue"]] <- "placenta"
  Middle3[["cohort"]] <- "Cohort7"
  Middle3[["group"]] <- "control"
  Middle3[["pregnancy_week"]] <- 18.2
  Middle3[["position"]] <- "VC"
  Middle4[["donor"]] <- "C7-P4"
  Middle4[["tissue"]] <- "placenta"
  Middle4[["cohort"]] <- "Cohort7"
  Middle4[["group"]] <- "control"
  Middle4[["pregnancy_week"]] <- 18.2
  Middle4[["position"]] <- "SC"
  Middle5[["donor"]] <- "C7-P5"
  Middle5[["tissue"]] <- "placenta"
  Middle5[["cohort"]] <- "Cohort7"
  Middle5[["group"]] <- "control"
  Middle5[["pregnancy_week"]] <- 17.6
  Middle5[["position"]] <- "VC"
  Middle6[["donor"]] <- "C7-P6"
  Middle6[["tissue"]] <- "placenta"
  Middle6[["cohort"]] <- "Cohort7"
  Middle6[["group"]] <- "control"
  Middle6[["pregnancy_week"]] <- 17.6
  Middle6[["position"]] <- "SC"
  Middle7[["donor"]] <- "C7-P7"
  Middle7[["tissue"]] <- "placenta"
  Middle7[["cohort"]] <- "Cohort7"
  Middle7[["group"]] <- "control"
  Middle7[["pregnancy_week"]] <- 23
  Middle7[["position"]] <- "VC"
  Middle8[["donor"]] <- "C7-P8"
  Middle8[["tissue"]] <- "placenta"
  Middle8[["cohort"]] <- "Cohort7"
  Middle8[["group"]] <- "control"
  Middle8[["pregnancy_week"]] <- 23
  Middle8[["position"]] <- "SC"
}
C7 <- merge(Middle1,c(Middle2,Middle3,Middle4,Middle5,Middle6,Middle7,Middle8))
C7[["pregnancy"]] <- "Middle"
save(C7,file = "./C7.Rdata")

#####QC#####
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

#####pre-process#####
merge.list <- list(C0,C1,C2,C3,C4,C5,C6,C7,C8)
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
var.features <- unname(obj = unlist(x = lapply(X = 1:length(x = merge.list), 
                                               FUN = function(x) VariableFeatures(object = merge.list[[x]], 
                                                                                  assay = "RNA"))))
write.csv(var.features,file = "./new_res/var.features.csv",row.names = F)
features <- SelectIntegrationFeatures(object.list = merge.list)
#####merge#####
merge.anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = features)
MFI_merge <- IntegrateData(anchorset = merge.anchors)
DefaultAssay(MFI_merge) <- "integrated"
#####downstream#####
MFI_merge <- ScaleData(MFI_merge, verbose = FALSE)
MFI_merge <- RunPCA(MFI_merge, npcs = 30, verbose = FALSE)
MFI_merge <- RunUMAP(MFI_merge, reduction = "pca", dims = 1:30)
save(MFI_merge,file = "MFI_merge.Rdata")