# Turn RNA rds into matrix

library(Seurat)
# library(Signac)

obj <- readRDS("../data/ML499M1-S1.rds")
obj

obj_assay = "RNA"
obj_layer = "data"

expr <- GetAssayData(obj, assay = obj_assay, layer = obj_layer)
expr <- t(expr)

dim(expr)

write.table(
  expr,
  file = "../data/ML499M1-S1_expression.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
cat("done")

