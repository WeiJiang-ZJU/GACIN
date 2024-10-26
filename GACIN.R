library(readr)
library(reshape2)
library(dplyr)
library(stringr)
library(tcltk)
##load data
load("./new_res/edges_melt_merge.Rdata")#all crosstalk terms
sample_info <- read.csv("./new_res/sample_info_20230319.csv")#MFI sample information
sample_info <- sample_info[which(sample_info$sampleName!=""),]
edges_melt_placenta <- edges_melt[,c(6,which(colnames(edges_melt)%in%sample_info$sampleName[which(sample_info$tissue=="placenta")]))]
#normal
MFI_normal <- edges_melt_placenta[,c(1,which(colnames(edges_melt_placenta)%in%sample_info$sampleName[which(sample_info$group=="normal")]))]
MFI_normal[is.na(MFI_normal)] <- 0
preg_normal <- sample_info$pregnancy_norm[match(colnames(MFI_normal)[-1],sample_info$sampleName)]
MFI_normal <- MFI_normal[which(rowSums(MFI_normal[,-1])!=0),]
#PE
PE <- edges_melt_placenta[,c(1,which(colnames(edges_melt_placenta)%in%sample_info$sampleName[which(sample_info$group=="PE")]))]
PE[is.na(PE)] <- 0
preg_PE <- sample_info$pregnancy_norm[match(colnames(PE)[-1],sample_info$sampleName)]
PE <- PE[which(rowSums(PE[,-1])!=0),]
#RPL
filename <- list.files("./new_res/GSE174399_cellphone_res/")
edges_PL <- read_tsv("./new_res/GSE174399_cellphone_res/A1_out_counts/significant_means.txt")
edges_PL <- edges_PL[!duplicated(edges_PL$interacting_pair),]
PL_melt <- melt(edges_PL,id.vars = "interacting_pair",variable.name = "cluster_pair",value.name = "A1",
                measure.vars = colnames(edges_PL)[13:ncol(edges_PL)])
PL_melt$cell_from <- unlist(lapply(strsplit(as.character(PL_melt$cluster_pair),"|",fixed = T), function(X) X[1]))
PL_melt$cell_to <- unlist(lapply(strsplit(as.character(PL_melt$cluster_pair),"|",fixed = T), function(X) X[2]))
PL_melt$crosstalk_term <- paste0(PL_melt$cell_from,"-",PL_melt$interacting_pair,"-",PL_melt$cell_to)
for (i in filename) {
  tmp <- read_tsv(paste0("./new_res/GSE174399_cellphone_res/",i,"/significant_means.txt"))
  tmp <- tmp[!duplicated(tmp$interacting_pair),]
  tmp_melt <- melt(tmp,id.vars = "interacting_pair",variable.name = "cluster_pair",value.name = str_extract(i,".*(?=_out_counts)"),
                   measure.vars = colnames(tmp)[(which(colnames(tmp)=="rank")+1):ncol(tmp)])
  tmp_melt$cell_from <- unlist(lapply(strsplit(as.character(tmp_melt$cluster_pair),"|",fixed = T), function(X) X[1]))
  tmp_melt$cell_to <- unlist(lapply(strsplit(as.character(tmp_melt$cluster_pair),"|",fixed = T), function(X) X[2]))
  tmp_melt$crosstalk_term <- paste0(tmp_melt$cell_from,"-",tmp_melt$interacting_pair,"-",tmp_melt$cell_to)
  PL_melt <- PL_melt %>% full_join(tmp_melt,by = "crosstalk_term")
}
PL_melt <- PL_melt[,which(colnames(PL_melt)%in%c("crosstalk_term",str_extract(filename[-21],".*(?=_out_counts)")))]
PL_melt[is.na(PL_melt)] <- 0
GSE174399_meta <- read.csv("./data/GSE174399_meta.csv")
preg_PL <- GSE174399_meta$pregnancy[match(colnames(PL_melt)[-2],GSE174399_meta$title)]
PL_melt <- PL_melt[which(rowSums(PL_melt[,-2])!=0),]
#PAS
filename <- list.files("./new_res/PAS/")
edges_PAS <- read_tsv("./new_res/PAS/GroupBP_out_counts/significant_means.txt")
edges_PAS <- edges_PAS[!duplicated(edges_PAS$interacting_pair),]
PAS_melt <- melt(edges_PAS,id.vars = "interacting_pair",variable.name = "cluster_pair",value.name = "GroupBP",
                 measure.vars = colnames(edges_PAS)[13:ncol(edges_PAS)])
PAS_melt$cell_from <- unlist(lapply(strsplit(as.character(PAS_melt$cluster_pair),"|",fixed = T), function(X) X[1]))
PAS_melt$cell_to <- unlist(lapply(strsplit(as.character(PAS_melt$cluster_pair),"|",fixed = T), function(X) X[2]))
PAS_melt$crosstalk_term <- paste0(PAS_melt$cell_from,"-",PAS_melt$interacting_pair,"-",PAS_melt$cell_to)
for (i in filename) {
  tmp <- read_tsv(paste0("./new_res/PAS/",i,"/significant_means.txt"))
  tmp <- tmp[!duplicated(tmp$interacting_pair),]
  tmp_melt <- melt(tmp,id.vars = "interacting_pair",variable.name = "cluster_pair",value.name = str_extract(i,".*(?=_out_counts)"),
                   measure.vars = colnames(tmp)[(which(colnames(tmp)=="rank")+1):ncol(tmp)])
  tmp_melt$cell_from <- unlist(lapply(strsplit(as.character(tmp_melt$cluster_pair),"|",fixed = T), function(X) X[1]))
  tmp_melt$cell_to <- unlist(lapply(strsplit(as.character(tmp_melt$cluster_pair),"|",fixed = T), function(X) X[2]))
  tmp_melt$crosstalk_term <- paste0(tmp_melt$cell_from,"-",tmp_melt$interacting_pair,"-",tmp_melt$cell_to)
  PAS_melt <- PAS_melt %>% full_join(tmp_melt,by = "crosstalk_term")
}
PAS_melt <- PAS_melt[,which(colnames(PAS_melt)%in%c("crosstalk_term",str_extract(filename,".*(?=_out_counts)")))]
PAS_melt[is.na(PAS_melt)] <- 0
PAS_melt <- PAS_melt[which(rowSums(PAS_melt[,-2])!=0),]

##correlation
MFI_normal_cor <- apply(MFI_normal[,-1],1, function(x) cor(x,preg_normal))
MFI_normal_crosstalk_term_cor <- data.frame(crosstalk_term = MFI_normal$crosstalk_term,
                                            cor = MFI_normal_cor)
MFI_normal_crosstalk_term_cor <- MFI_normal_crosstalk_term_cor[order(MFI_normal_crosstalk_term_cor$cor,decreasing = T),]
write.csv(MFI_normal_crosstalk_term_cor,file = "./new_res/clock/all_merge/MFI_normal_crosstalk_term_cor(edit_placenta).csv",row.names = F)
MFI_normal_cor_test <- apply(MFI_normal[,-1],1, function(x) cor.test(x,preg_normal))
MFI_normal_cor_test <- lapply(MFI_normal_cor_test,FUN = function(x) x$p.value)
MFI_normal_cor_test <- unlist(MFI_normal_cor_test)
MFI_normal_crosstalk_term_cor_test <- data.frame(crosstalk_term = MFI_normal$crosstalk_term,
                                                 cor = MFI_normal_cor,
                                                 pval = MFI_normal_cor_test)
MFI_normal_crosstalk_term_cor_test <- MFI_normal_crosstalk_term_cor_test[order(MFI_normal_crosstalk_term_cor_test$cor,decreasing = T),]

##EPL-DE
MFI_normal_early <- MFI_normal[,which(colnames(MFI_normal)%in%c("crosstalk_term",sample_info$sampleName[which(sample_info$time=="Early")]))]
EPL <- PL_melt[,c(2,1,3:11)]
EPL <- EPL[which(rowSums(EPL[,-1])!=0),]
RPL <- PL_melt[,c(2,12:16)]
RPL <- RPL[which(rowSums(RPL[,-1])!=0),]
PL_normal <- PL_melt[,c(2,17:21)]
PL_normal <- PL_normal[which(rowSums(PL_normal[,-1])!=0),]
EPL_normal_mix <- EPL %>% full_join(MFI_normal_early,by = "crosstalk_term")
EPL_normal_mix[is.na(EPL_normal_mix)] <- 0
shap_test_EPL <- list()
for(i in 1:nrow(EPL_normal_mix[,2:11])){
  if(sum(EPL_normal_mix[i,2:11])==0){shap_test_EPL[[i]] <- 1}
  else{
    shap_test_EPL[[i]] <- shapiro.test(as.numeric(EPL_normal_mix[i,2:11]))
    shap_test_EPL[[i]] <- shap_test_EPL[[i]]$p.value
  }
}
shap_test_EPL_if <- unlist(shap_test_EPL)
shap_test_norm_early <- list()
for(i in 1:nrow(EPL_normal_mix[,12:33])){
  if(sum(EPL_normal_mix[i,12:33])==0){shap_test_norm_early[[i]] <- 1}
  else{
    shap_test_norm_early[[i]] <- shapiro.test(as.numeric(EPL_normal_mix[i,12:33]))
    shap_test_norm_early[[i]] <- shap_test_norm_early[[i]]$p.value
  }
}
shap_test_norm_early_if <- unlist(shap_test_norm_early)
EPL_test_res <- list()
shap_test_norm_if <- shap_test_norm_early_if >= 0.05
shap_test_EPL_if <- shap_test_EPL_if >= 0.05
indices <- which(shap_test_norm_if & shap_test_EPL_if)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- t.test(as.numeric(EPL_normal_mix[indices[i],2:11]),as.numeric(EPL_normal_mix[indices[i],12:33]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
EPL_test_res <- data.frame(pval = pvals, method = methods)
row.names(EPL_test_res) <- EPL_normal_mix$crosstalk_term[indices]
indices <- which((shap_test_norm_if & shap_test_EPL_if)==F)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- wilcox.test(as.numeric(EPL_normal_mix[indices[i],2:11]),as.numeric(EPL_normal_mix[indices[i],12:33]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
EPL_test_res_wilcox <- data.frame(pval = pvals, method = methods)
row.names(EPL_test_res_wilcox) <- EPL_normal_mix$crosstalk_term[indices]
EPL_test_res <- rbind(EPL_test_res,EPL_test_res_wilcox)

##RPL-DE
MFI_normal_early <- MFI_normal[,which(colnames(MFI_normal)%in%c("crosstalk_term",sample_info$sampleName[which(sample_info$time=="Early")]))]
RPL <- PL_melt[,c(2,12:16)]
RPL <- RPL[which(rowSums(RPL[,-1])!=0),]
RPL_normal_mix <- RPL %>% full_join(MFI_normal_early,by = "crosstalk_term")
RPL_normal_mix[is.na(RPL_normal_mix)] <- 0
shap_test_RPL <- list()
for(i in 1:nrow(RPL_normal_mix[,2:6])){
  if(sum(RPL_normal_mix[i,2:6])==0){shap_test_RPL[[i]] <- 1}
  else{
    shap_test_RPL[[i]] <- shapiro.test(as.numeric(RPL_normal_mix[i,2:6]))
    shap_test_RPL[[i]] <- shap_test_RPL[[i]]$p.value
  }
}
shap_test_RPL_if <- unlist(shap_test_RPL)
shap_test_norm_early <- list()
for(i in 1:nrow(RPL_normal_mix[,7:28])){
  if(sum(RPL_normal_mix[i,7:28])==0){shap_test_norm_early[[i]] <- 1}
  else{
    shap_test_norm_early[[i]] <- shapiro.test(as.numeric(RPL_normal_mix[i,7:28]))
    shap_test_norm_early[[i]] <- shap_test_norm_early[[i]]$p.value
  }
}
shap_test_norm_early_if <- unlist(shap_test_norm_early)
shap_test_norm_if <- shap_test_norm_early_if >= 0.05
shap_test_RPL_if <- shap_test_RPL_if >= 0.05
indices <- which(shap_test_norm_if & shap_test_RPL_if)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- t.test(as.numeric(RPL_normal_mix[indices[i],2:6]),as.numeric(RPL_normal_mix[indices[i],7:28]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
RPL_test_res <- data.frame(pval = pvals, method = methods)
row.names(RPL_test_res) <- RPL_normal_mix$crosstalk_term[indices]
indices <- which((shap_test_norm_if & shap_test_RPL_if)==F)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- wilcox.test(as.numeric(RPL_normal_mix[indices[i],2:6]),as.numeric(RPL_normal_mix[indices[i],7:28]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
RPL_test_res_wilcox <- data.frame(pval = pvals, method = methods)
row.names(RPL_test_res_wilcox) <- RPL_normal_mix$crosstalk_term[indices]
RPL_test_res <- rbind(RPL_test_res,RPL_test_res_wilcox)

##PAS-DE
MFI_normal_middle <- MFI_normal[,which(colnames(MFI_normal)%in%c("crosstalk_term",sample_info$sampleName[which(sample_info$time=="Middle")]))]
MFI_normal_middle <- MFI_normal_middle[which(rowSums(MFI_normal_middle[,-1])!=0),]
PAS <- PAS_melt[,-4]
PAS_normal_mix <- PAS %>% full_join(MFI_normal_middle,by = "crosstalk_term")
PAS_normal_mix[is.na(PAS_normal_mix)] <- 0
PAS_normal_mix <- PAS_normal_mix[,c(2,1,3:20)]
shap_test_PAS <- list()
for(i in 1:nrow(PAS_normal_mix[,2:4])){
  if(sum(PAS_normal_mix[i,2:4])==0){shap_test_PAS[[i]] <- 1}
  else{
    shap_test_PAS[[i]] <- shapiro.test(as.numeric(PAS_normal_mix[i,2:4]))
    shap_test_PAS[[i]] <- shap_test_PAS[[i]]$p.value
  }
}
shap_test_PAS_if <- unlist(shap_test_PAS)
shap_test_norm_middle <- list()
for(i in 1:nrow(PAS_normal_mix[,5:20])){
  if(sum(PAS_normal_mix[i,5:20])==0){shap_test_norm_middle[[i]] <- 1}
  else{
    shap_test_norm_middle[[i]] <- shapiro.test(as.numeric(PAS_normal_mix[i,5:20]))
    shap_test_norm_middle[[i]] <- shap_test_norm_middle[[i]]$p.value
  }
}
shap_test_norm_middle_if <- unlist(shap_test_norm_middle)
shap_test_norm_if <- shap_test_norm_middle_if >= 0.05
shap_test_PAS_if <- shap_test_PAS_if >= 0.05
indices <- which(shap_test_norm_if & shap_test_PAS_if)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- t.test(as.numeric(PAS_normal_mix[indices[i],2:4]),as.numeric(PAS_normal_mix[indices[i],5:20]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
PAS_test_res <- data.frame(pval = pvals, method = methods)
row.names(PAS_test_res) <- PAS_normal_mix$crosstalk_term[indices]
indices <- which((shap_test_norm_if & shap_test_PAS_if)==F)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- wilcox.test(as.numeric(PAS_normal_mix[indices[i],2:4]),as.numeric(PAS_normal_mix[indices[i],5:20]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
PAS_test_res_wilcox <- data.frame(pval = pvals, method = methods)
row.names(PAS_test_res_wilcox) <- PAS_normal_mix$crosstalk_term[indices]
PAS_test_res <- rbind(PAS_test_res,PAS_test_res_wilcox)

##PE-DE
MFI_normal_late <- MFI_normal[,which(colnames(MFI_normal)%in%c("crosstalk_term",sample_info$sampleName[which(sample_info$time=="Late")]))]
MFI_normal_late <- MFI_normal_late[which(rowSums(MFI_normal_late[,-1])!=0),]
PE_normal_mix <- PE %>% full_join(MFI_normal_late,by = "crosstalk_term")
PE_normal_mix[is.na(PE_normal_mix)] <- 0
shap_test_PE <- list()
for(i in 1:nrow(PE_normal_mix[,2:9])){
  if(sum(PE_normal_mix[i,2:9])==0){shap_test_PE[[i]] <- 1}
  else{
    shap_test_PE[[i]] <- shapiro.test(as.numeric(PE_normal_mix[i,2:9]))
    shap_test_PE[[i]] <- shap_test_PE[[i]]$p.value
  }
}
shap_test_PE_if <- unlist(shap_test_PE)
shap_test_norm_late <- list()
for(i in 1:nrow(PE_normal_mix[,10:19])){
  if(sum(PE_normal_mix[i,10:19])==0){shap_test_norm_late[[i]] <- 1}
  else{
    shap_test_norm_late[[i]] <- shapiro.test(as.numeric(PE_normal_mix[i,10:19]))
    shap_test_norm_late[[i]] <- shap_test_norm_late[[i]]$p.value
  }
}
shap_test_norm_late_if <- unlist(shap_test_norm_late)
shap_test_norm_if <- shap_test_norm_late_if >= 0.05
shap_test_PE_if <- shap_test_PE_if >= 0.05
indices <- which(shap_test_norm_if & shap_test_PE_if)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- t.test(as.numeric(PE_normal_mix[indices[i],2:9]),as.numeric(PE_normal_mix[indices[i],10:19]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
PE_test_res <- data.frame(pval = pvals, method = methods)
row.names(PE_test_res) <- PE_normal_mix$crosstalk_term[indices]
indices <- which((shap_test_norm_if & shap_test_PE_if)==F)
pvals <- numeric(length(indices))
methods <- character(length(indices))
for(i in seq_along(indices)){
  t <- wilcox.test(as.numeric(PE_normal_mix[indices[i],2:9]),as.numeric(PE_normal_mix[indices[i],10:19]),paired = F)
  pvals[i] <- t$p.value
  methods[i] <- t$method
}
PE_test_res_wilcox <- data.frame(pval = pvals, method = methods)
row.names(PE_test_res_wilcox) <- PE_normal_mix$crosstalk_term[indices]
PE_test_res <- rbind(PE_test_res,PE_test_res_wilcox)

#####GACIN-----
normal_sig <- MFI_normal_crosstalk_term_cor_test[which(abs(MFI_normal_crosstalk_term_cor_test$pval)<0.05),]
crosstalk_term_normal_top20 <- normal_sig[c(1:10,15634:15643),]
mtx_normal <- crosstalk_term_normal_top20 %>% left_join(MFI_normal,by = "crosstalk_term")
mtx_normal <- mtx_normal[,-c(2:3)]
mtx_normal[is.na(mtx_normal)] <- 0
fit_mtx_top <- as.data.frame(t(mtx_normal[,-1]))
colnames(fit_mtx_top) <- mtx_normal$crosstalk_term
fit_mtx_top$pregnancy <- preg_normal
fit_top <- lm(formula = pregnancy ~ ., 
              data = fit_mtx_top)
summary(fit_top)
library(gvlma)
gvmodel <- gvlma(fit_top)
summary(gvmodel)
library(MASS)
stepAIC(fit_top,direction = "backward")
fit_backward_reg <- lm(pregnancy ~ `SCT-CD55_ADGRE5-dM2` + `dM1-ICAM3_CD209-HB` + 
                         `dM3-CD74_MIF-CD8+ naive T` + `dM2-LAIR1_LILRB4-HB` + `dNK3-CD94:NKG2A_HLA-E-CD8+ naive T` + 
                         `CD8+ naive T-HLA-E_KLRC1-dNK3` + `fFB1-MDK_LRP1-VCT` + `fFB1-MDK_LRP1-fFB1` + 
                         `fFB1-MDK_LRP1-dM3`,data = fit_mtx_top)
summary(fit_backward_reg)
gvmodel <- gvlma(fit_backward_reg)
summary(gvmodel)

library(leaps)
fit_leaps_scale <- regsubsets(pregnancy ~ .,data=fit_mtx_top,nbest = 100,really.big = T)
plot(fit_leaps_scale,scale = "adjr2")
colnames(fit_mtx_top)[c(5,6,7,10,18,19)]
fit_leaps_reg_top_cor <- lm(formula = pregnancy ~ 
                              `dM3-CD74_MIF-CD8+ naive T` +
                              `dM2-LAIR1_LILRB4-HB` + `DC2-ICAM3_CD209-HB` +
                              `dM3-GRN_SORT1-Endo (m)`+
                              `fFB1-MDK_LRP1-VCT`+
                              `fFB1-MDK_LRP1-dM3`,
                            data = fit_mtx_top)
summary(fit_leaps_reg_top_cor)
gvmodel <- gvlma(fit_leaps_reg_top_cor)
summary(gvmodel)
###k-fold cross validation
library(caret)
set.seed(121)
custom_summary_function <- function(data, lev = NULL, model = NULL) {
  predictions <- data$pred
  return(predictions)
}
train.control <- trainControl(method = "repeatedcv",
                              number = 6, repeats = 200)
model_fold <- train(pregnancy ~ 
                      `dM3-CD74_MIF-CD8+ naive T` +
                      `dM2-LAIR1_LILRB4-HB` + `DC2-ICAM3_CD209-HB` +
                      `dM3-GRN_SORT1-Endo (m)`+
                      `fFB1-MDK_LRP1-VCT`+
                      `fFB1-MDK_LRP1-dM3`,
                    data = fit_mtx_top, method = "lm",
                    trControl = train.control)
print(model_fold)
###visualization
library(ggstatsplot)
library(ggthemes)
ggcoefstats(model_fold$finalModel,stats.label.color = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:7],
            point.args = list(size = 3, color = "#00409a", na.rm = TRUE),
            errorbar.args = list(height = 0, color = "#00409a", na.rm = TRUE),
            vline.args = list(linewidth = 0.5, linetype = "dashed")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank())