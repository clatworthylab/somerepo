```R
library(miloR)
library(igraph)
library(BiocParallel)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
setwd('/lustre/scratch117/cellgen/team297/kt16/Ziad/scanpy')

load('h5ad/R/classical_mono_milo_prep_1.5.RData')
# 1.5 MIU vs Untreated (B)
mylox1 <- Milo(adata_no_knnx1)
milo_graphx1 <- buildFromAdjacency(knn_adjacencyx1, k=50, is.binary=TRUE)
graph(mylox1) <- graph(milo_graphx1)
mylox1 <- buildGraph(mylox1, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylox1 <- makeNhoods(mylox1, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylox1 <- countCells_(mylox1, meta.data = data.frame(colData(mylox1)), samples="sample_id")
mylox1 <- calcNhoodDistance(mylox1, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylox1 <- buildNhoodGraph(mylox1)
saveRDS(mylox1, 'h5ad/R/milo_classical_mono_1.5.RDS')

load('h5ad/R/classical_mono_milo_prep_2.5.RData')
# 2.5 MIU vs Untreated (B)
mylox2 <- Milo(adata_no_knnx2)
milo_graphx2 <- buildFromAdjacency(knn_adjacencyx2, k=50, is.binary=TRUE)
graph(mylox2) <- graph(milo_graphx2)
mylox2 <- buildGraph(mylox2, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylox2 <- makeNhoods(mylox2, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylox2 <- countCells_(mylox2, meta.data = data.frame(colData(mylox2)), samples="sample_id")
mylox2 <- calcNhoodDistance(mylox2, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylox2 <- buildNhoodGraph(mylox2)
saveRDS(mylox2, 'h5ad/R/milo_classical_mono_2.5.RDS')

```
### so people in teichlab are concerned that there may be too many covariates so let's remove lymph perhaps since this is the one that probably the most alarming.
```R
library(miloR)
library(igraph)
library(BiocParallel)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(lmerTest)
library(pbmcapply)
setwd('/lustre/scratch117/cellgen/team297/kt16/Ziad/scanpy')

metadata = read.table('../sampleinfo.txt', header = 1)
metadata$timepoint = factor(metadata$timepoint, level= c('pre', 'post'))
metadata$treatment = factor(metadata$treatment, level= c('Saline', '1.5MIU', '2.5MIU'))
metadata$treatment_group_1 = factor(metadata$treatment_group_1, level= c('untreated', '1.5MIU', '2.5MIU'))
metadata$treatment_group_1_ordered = ordered(metadata$treatment_group_1, level= c('untreated', '1.5MIU', '2.5MIU'))
metadata$treatment_group_2 = factor(metadata$treatment_group_2, level= c('untreated', 'treated'))
# scale otherwise it returns lots of warnings
metadata$peak_trop <- scale(metadata$peak_trop)
metadata$age <- scale(metadata$age)
metadata$Lymph <- scale(metadata$Lymph)

metadata1 <- metadata %>% filter(treatment %in% c('Saline', '1.5MIU'))
metadata1$treatment <- droplevels(metadata1$treatment)

metadata2 <- metadata %>% filter(treatment %in% c('Saline', '2.5MIU'))
metadata2$treatment <- droplevels(metadata2$treatment)

mylox1 <- readRDS('h5ad/R/milo_classical_mono_1.5.RDS')
mylox2 <- readRDS('h5ad/R/milo_classical_mono_2.5.RDS')

resx1b <- as.list(1:nrow(nhoodCounts(mylox1)))
resx1b <- pbmclapply(resx1b, function(x) glmer.nb(nhoodCounts(mylox1)[x,] ~ sex + age + peak_trop + (1|study_id) + treatment_group_1, data = metadata1, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())
res_summaryx1b <- pbmclapply(resx1b, summary)
save(res_summaryx1b, resx1b, metadata1, file = 'h5ad/R/model_1/milo_classical_mono_1.5_results.RData')

resx2b <- as.list(1:nrow(nhoodCounts(mylox2)))
resx2b <- pbmclapply(resx2b, function(x) glmer.nb(nhoodCounts(mylox2)[x,] ~ sex + age + peak_trop + (1|study_id) + treatment_group_1, data = metadata2, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())
res_summaryx2b <- pbmclapply(resx2b, summary)
save(res_summaryx2b, resx2b, metadata2, file = 'h5ad/R/model_1/milo_classical_mono_2.5_results.RData')

```

# Export results
```R
library(miloR)
library(pbmcapply)
setwd('/lustre/scratch117/cellgen/team297/kt16/Ziad/scanpy')

mylox1 <- readRDS('h5ad/R/milo_classical_mono_1.5.RDS')
mylox2 <- readRDS('h5ad/R/milo_classical_mono_2.5.RDS')

load('h5ad/R/model_2/milo_classical_mono_1.5_results.RData')
load('h5ad/R/model_2/milo_classical_mono_2.5_results.RData')

b1.5_pvalue <- pbmclapply(res_summaryx1b, function(x) x$coefficients[5,4])
b2.5_pvalue <- pbmclapply(res_summaryx2b, function(x) x$coefficients[5,4])

b1.5_beta <- pbmclapply(res_summaryx1b, function(x) x$coefficients[5,1])
b2.5_beta <- pbmclapply(res_summaryx2b, function(x) x$coefficients[5,1])

da.res_b1.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylox1))), beta = do.call(c, b1.5_beta), PValue = do.call(c, b1.5_pvalue))
da.res_b2.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylox2))), beta = do.call(c, b2.5_beta), PValue = do.call(c, b2.5_pvalue))

# da.res$Nhood <- as.numeric(rownames(da.res))
message("Computing SpatialFDR")
da.res_b1.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylox1), graph=miloR::graph(mylox1),
	weighting='k-distance', pvalues=da.res_b1.5$PValue,
	indices=nhoodIndex(mylox1), distances=nhoodDistances(mylox1),
	reduced.dimensions=reducedDim(mylox1, 'X_pca_harmony'), k = 50)
da.res_b1.5$Diff <- sign(da.res_b1.5$beta)
da.res_b1.5$Diff[da.res_b1.5$SpatialFDR >= 0.1] <- 0

da.res_b2.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylox2), graph=miloR::graph(mylox2),
	weighting='k-distance', pvalues=da.res_b2.5$PValue,
	indices=nhoodIndex(mylox2), distances=nhoodDistances(mylox2),
	reduced.dimensions=reducedDim(mylox2, 'X_pca_harmony'), k = 50)
da.res_b2.5$Diff <- sign(da.res_b2.5$beta)
da.res_b2.5$Diff[da.res_b2.5$SpatialFDR >= 0.1] <- 0

saveRDS(da.res_b1.5, file = 'h5ad/R/model_1/milo_results_1.5vsuntreated_classical_mono.RDS')
saveRDS(da.res_b2.5, file = 'h5ad/R/model_1/milo_results_2.5vsuntreated_classical_mono.RDS')

```
