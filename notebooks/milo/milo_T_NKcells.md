```R
library(miloR)
library(igraph)
library(BiocParallel)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
setwd('/lustre/scratch117/cellgen/team297/kt16/Ziad/scanpy')

load('h5ad/R/milo_prep_1.5.RData')
# 1.5 MIU vs Untreated (all)
mylo1 <- Milo(adata_no_knn1)
milo_graph1 <- buildFromAdjacency(knn_adjacency1, k=50, is.binary=TRUE)
graph(mylo1) <- graph(milo_graph1)
mylo1 <- buildGraph(mylo1, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylo1 <- makeNhoods(mylo1, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylo1 <- countCells_(mylo1, meta.data = data.frame(colData(mylo1)), samples="sampleid")
mylo1 <- calcNhoodDistance(mylo1, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylo1 <- buildNhoodGraph(mylo1)
saveRDS(mylo1, 'h5ad/R/milo_full_1.5.RDS')

# 1.5 MIU vs Untreated (TNK)
mylo1 <- Milo(adata_no_knnx1)
milo_graphx1 <- buildFromAdjacency(knn_adjacencyx1, k=50, is.binary=TRUE)
graph(mylo1) <- graph(milo_graphx1)
mylo1 <- buildGraph(mylo1, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylo1 <- makeNhoods(mylo1, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylo1 <- countCells_(mylo1, meta.data = data.frame(colData(mylo1)), samples="sampleid")
mylo1 <- calcNhoodDistance(mylo1, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylo1 <- buildNhoodGraph(mylo1)
saveRDS(mylo1, 'h5ad/R/milo_TNK_1.5.RDS')


load('h5ad/R/milo_prep_2.5.RData')
# 2.5 MIU vs Untreated (all)
mylo2 <- Milo(adata_no_knn2)
milo_graph2 <- buildFromAdjacency(knn_adjacency2, k=50, is.binary=TRUE)
graph(mylo2) <- graph(milo_graph2)
mylo2 <- buildGraph(mylo2, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylo2 <- makeNhoods(mylo2, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylo2 <- countCells_(mylo2, meta.data = data.frame(colData(mylo2)), samples="sampleid")
mylo2 <- calcNhoodDistance(mylo2, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylo2 <- buildNhoodGraph(mylo2)
saveRDS(mylo2, 'h5ad/R/milo_full_2.5.RDS')

# 2.5 MIU vs Untreated (TNK)
mylo2 <- Milo(adata_no_knnx2)
milo_graphx2 <- buildFromAdjacency(knn_adjacencyx2, k=50, is.binary=TRUE)
graph(mylo2) <- graph(milo_graphx2)
mylo2 <- buildGraph(mylo2, k=50, d=50, reduced.dim="X_pca_harmony", BPPARAM = MulticoreParam(progressbar = TRUE))
mylo2 <- makeNhoods(mylo2, prop = 0.1, k = 50, d=50, reduced_dim="X_pca_harmony")
mylo2 <- countCells_(mylo2, meta.data = data.frame(colData(mylo2)), samples="sampleid")
mylo2 <- calcNhoodDistance(mylo2, d=50, reduced.dim = 'X_pca_harmony', use.assay = 'counts')
mylo2 <- buildNhoodGraph(mylo2)
saveRDS(mylo2, 'h5ad/R/milo_TNK_2.5.RDS')

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

mylo1 <- readRDS('h5ad/R/milo_full_1.5.RDS')
mylo2 <- readRDS('h5ad/R/milo_full_2.5.RDS')

res1 <- as.list(1:nrow(nhoodCounts(mylo1)))
res1 <- pbmclapply(res1, function(x) glmer.nb(nhoodCounts(mylo1)[x,] ~ sex + age + peak_trop + (1|study_id) + offset(log(colSums(nhoodCounts(mylo1)))) + treatment_group_1, data = metadata1, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())

res_summary1 <- pbmclapply(res1, summary)
save(res_summary1, res1, metadata1, file = 'h5ad/R/model_1/milo_full_1.5_results.RData')

res2 <- as.list(1:nrow(nhoodCounts(mylo2)))
res2 <- pbmclapply(res2, function(x) glmer.nb(nhoodCounts(mylo2)[x,] ~ sex + age + peak_trop + (1|study_id) + offset(log(colSums(nhoodCounts(mylo2)))) + treatment_group_1, data = metadata2, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())

res_summary2 <- pbmclapply(res2, summary)
save(res_summary2, res2, metadata2, file = 'h5ad/R/model_1/milo_full_2.5_results.RData')

## just T and NK cells
mylo1 <- readRDS('h5ad/R/milo_TNK_1.5.RDS')
mylo2 <- readRDS('h5ad/R/milo_TNK_2.5.RDS')

res1 <- as.list(1:nrow(nhoodCounts(mylo1)))
res1 <- pbmclapply(res1, function(x) glmer.nb(nhoodCounts(mylo1)[x,] ~ sex + age + peak_trop + (1|study_id) + offset(log(colSums(nhoodCounts(mylo1)))) + treatment_group_1, data = metadata1, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())
res_summary1 <- pbmclapply(res1, summary)
save(res_summary1, res1, metadata1, file = 'h5ad/R/model_1/milo_TNK_1.5_results.RData')

res2 <- as.list(1:nrow(nhoodCounts(mylo2)))
res2 <- pbmclapply(res2, function(x) glmer.nb(nhoodCounts(mylo2)[x,] ~ sex + age + peak_trop + (1|study_id) + offset(log(colSums(nhoodCounts(mylo2)))) + treatment_group_1, data = metadata2, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-2), optimizer="bobyqa", tol=1e-02)), mc.cores = parallel::detectCores())
res_summary2 <- pbmclapply(res2, summary)
save(res_summary2, res2, metadata2, file = 'h5ad/R/model_1/milo_TNK_2.5_results.RData')

```

```R
library(miloR)
library(pbmcapply)
setwd('/lustre/scratch117/cellgen/team297/kt16/Ziad/scanpy')

mylo1 <- readRDS('h5ad/R/milo_TNK_1.5.RDS')
mylo2 <- readRDS('h5ad/R/milo_TNK_2.5.RDS')
load('h5ad/R/model_1/milo_TNK_1.5_results.RData')
load('h5ad/R/model_1/milo_TNK_2.5_results.RData')
tnk1.5_pvalue <- pbmclapply(res_summary1, function(x) x$coefficients[5,4])
tnk2.5_pvalue <- pbmclapply(res_summary2, function(x) x$coefficients[5,4])
tnk1.5_beta <- pbmclapply(res_summary1, function(x) x$coefficients[5,1])
tnk2.5_beta <- pbmclapply(res_summary2, function(x) x$coefficients[5,1])
da.res_tnk1.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylo1))), beta = do.call(c, tnk1.5_beta), PValue = do.call(c, tnk1.5_pvalue))
da.res_tnk2.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylo2))), beta = do.call(c, tnk2.5_beta), PValue = do.call(c, tnk2.5_pvalue))
message("Computing SpatialFDR")
da.res_tnk1.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylo1), graph=miloR::graph(mylo1),
	weighting='k-distance', pvalues=da.res_tnk1.5$PValue,
	indices=nhoodIndex(mylo1), distances=nhoodDistances(mylo1),
	reduced.dimensions=reducedDim(mylo1, 'X_pca_harmony'), k = 50)
da.res_tnk1.5$Diff <- sign(da.res_tnk1.5$beta)
da.res_tnk1.5$Diff[da.res_tnk1.5$SpatialFDR >= 0.1] <- 0
da.res_tnk2.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylo2), graph=miloR::graph(mylo2),
	weighting='k-distance', pvalues=da.res_tnk2.5$PValue,
	indices=nhoodIndex(mylo2), distances=nhoodDistances(mylo2),
	reduced.dimensions=reducedDim(mylo2, 'X_pca_harmony'), k = 50)
da.res_tnk2.5$Diff <- sign(da.res_tnk2.5$beta)
da.res_tnk2.5$Diff[da.res_tnk2.5$SpatialFDR >= 0.1] <- 0
saveRDS(da.res_tnk1.5, file = 'h5ad/R/model_1/milo_results_1.5vsuntreated_TNK.RDS')
saveRDS(da.res_tnk2.5, file = 'h5ad/R/model_1/milo_results_2.5vsuntreated_TNK.RDS')

mylo1 <- readRDS('h5ad/R/milo_full_1.5.RDS')
mylo2 <- readRDS('h5ad/R/milo_full_2.5.RDS')
load('h5ad/R/model_1/milo_full_1.5_results.RData')
load('h5ad/R/model_1/milo_full_2.5_results.RData')
full1.5_pvalue <- pbmclapply(res_summary1, function(x) x$coefficients[5,4])
full2.5_pvalue <- pbmclapply(res_summary2, function(x) x$coefficients[5,4])
full1.5_beta <- pbmclapply(res_summary1, function(x) x$coefficients[5,1])
full2.5_beta <- pbmclapply(res_summary2, function(x) x$coefficients[5,1])
da.res_full1.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylo1))), beta = do.call(c, full1.5_beta), PValue = do.call(c, full1.5_pvalue))
da.res_full2.5 = data.frame(Nhood = as.numeric(1:nrow(nhoodCounts(mylo2))), beta = do.call(c, full2.5_beta), PValue = do.call(c, full2.5_pvalue))
da.res_full1.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylo1), graph=miloR::graph(mylo1),
	weighting='k-distance', pvalues=da.res_full1.5$PValue,
	indices=nhoodIndex(mylo1), distances=nhoodDistances(mylo1),
	reduced.dimensions=reducedDim(mylo1, 'X_pca_harmony'), k = 50)
da.res_full1.5$Diff <- sign(da.res_full1.5$beta)
da.res_full1.5$Diff[da.res_full1.5$SpatialFDR >= 0.1] <- 0
da.res_full2.5$SpatialFDR <- graphSpatialFDR(x.nhoods=nhoods(mylo2), graph=miloR::graph(mylo2),
	weighting='k-distance', pvalues=da.res_full2.5$PValue,
	indices=nhoodIndex(mylo2), distances=nhoodDistances(mylo2),
	reduced.dimensions=reducedDim(mylo2, 'X_pca_harmony'), k = 50)
da.res_full2.5$Diff <- sign(da.res_full2.5$beta)
da.res_full2.5$Diff[da.res_full2.5$SpatialFDR >= 0.1] <- 0
saveRDS(da.res_full1.5, file = 'h5ad/R/model_1/milo_results_1.5vsuntreated_full.RDS')
saveRDS(da.res_full2.5, file = 'h5ad/R/model_1/milo_results_2.5vsuntreated_full.RDS')
```