# Description: Extract highly variable genes of zebrafish embryonic (ZB) scRNA-seq data.
# Author: Jiaqi Zhang <jiaqi_zhang2@brown.edu>
# Reference:
#   https://www.science.org/doi/10.1126/science.aar3131
#   https://github.com/farrellja/URD
#   https://singlecell.broadinstitute.org/single_cell/study/SCP162/
library(Seurat)
library(stringr)

# -----
# Load data
print("Loading zebrafish data...")
data_obj <- readRDS("URD_Zebrafish_Object.rds")
data_count <- data_obj@count.data # gene x cell
meta_df <- data_obj@meta
print(sprintf("Data shape (gene x cell): %d x %d", dim(data_count)[1], dim(data_count)[2]))
# -----
# Split data by traing and testing sets
split_type <- "three_interpolation" # two_forecasting, three_interpolation, remove_recovery
cell_tp <- meta_df$stage.nice
unique_tp <- unique(cell_tp)

train_tps <- c(1, 2, 3, 4, 6, 8, 10, 11, 12)
test_tps <- c(5, 7, 9)

print(sprintf("Num tps: %d", length(unique_tp)))
print("Train tps:")
print(train_tps)
print("Test tps:")
print(test_tps)
train_data <- data_count[, which(cell_tp %in% unique_tp[train_tps])]
# test_data <- data_count[, which(cell_tp %in% unique_tp[test_tps])]
print(sprintf("Train data shape (gene x cell): %d x %d", dim(train_data)[1], dim(train_data)[2]))
# print(sprintf("Test data  shape (gene x cell): %d x %d", dim(test_data)[1], dim(test_data)[2]))
# -----
# Select highly variables based on training data
train_obj <- CreateSeuratObject(counts=train_data)
train_obj <- FindVariableFeatures(NormalizeData(train_obj), selection.method = "vst", nfeatures = 2000) # use log-normalized data
hvgs <- VariableFeatures(train_obj)
hvg_diff <- hvgs[!(hvgs %in% row.names(data_count))]
if (length(hvg_diff) > 0){ # solve the problem that replacing "_" with "-"
  print("Tune hvg list...")
  hvgs <- c(hvgs[hvgs %in% row.names(data_count)], str_replace(hvg_diff, "-", "_"))
}
data_count_hvg <- data_count[hvgs, ]
print(sprintf("HVG data shape (gene x cell): %d x %d", dim(data_count_hvg)[1], dim(data_count_hvg)[2]))

write.csv(as.matrix(t(data_count_hvg)), "count_data-hvg.csv")
write.csv(hvgs, "var_genes_list.csv")
write.csv(meta_df, "meta_data.csv")
