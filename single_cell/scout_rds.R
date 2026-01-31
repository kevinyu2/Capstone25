library(Seurat)
library(Signac)

obj <- readRDS("CESC_CE348E1-S1K1.rds")
class(obj)
obj


obj_assay = "peaks"
obj_layer = "counts"

counts <- GetAssayData(
  object = obj,
  assay = obj_assay,
  layer = obj_layer
)

rownames(counts) # should be genes or ATAC locations
colnames(counts) # should be cells

cell_test_name <- colnames(counts)[1]
cell_test <- counts[, cell_test_name]

value_counts <- table(cell_test)
print(value_counts)









# Extract the numbers before the `#` (sample IDs)
sample_ids <- substr(colnames(counts), 1, 3)


# Get all unique sample IDs
unique_sample_ids <- unique(sample_ids)

# Print the unique sample IDs
print(unique_sample_ids)