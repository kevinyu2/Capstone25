# Turn RNA rds into matrix

library(Seurat)
# library(Signac)

obj <- readRDS("./data/HT268B1-Th1H3.rds")
obj


obj_assay = "SCT"
obj_layer = "scale.data"

expr <- GetAssayData(obj, assay = obj_assay, layer = obj_layer)
expr <- t(expr)

dim(expr)

write.table(
  expr,
  file = "./data/sct_scale_HT268B1-Th1H3.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
cat("done")

obj
obj_assay = "RNA"
obj_layer = "counts"

expr <- GetAssayData(obj, assay = obj_assay, layer = obj_layer)
expr <- t(expr)

dim(expr)

write.table(
  expr,
  file = "./data/expression_counts_HT268B1-Th1H3.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
cat("done")

