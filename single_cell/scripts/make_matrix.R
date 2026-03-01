# Turn RNA rds into matrix

library(Seurat)
# library(Signac)

obj <- readRDS("./data/ML499M1-S1.rds")
obj

obj.

obj_assay = "SCT"
obj_layer = "scale.data"

expr <- GetAssayData(obj, assay = obj_assay, layer = obj_layer)
expr <- t(expr)

dim(expr)

write.table(
  expr,
  file = "./data/sct_scale_ML499M1-S1.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
cat("done")

