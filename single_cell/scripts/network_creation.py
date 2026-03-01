import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import cosine_similarity

# TODO: scaling looks bad, but maybe inter-percentile range? That would essentially be like 95% - 5% as the variance, this might capture bimodality and ecDNA better


################################################
## Settings
################################################
cell_by_gene = pd.read_csv('../data/count_matrix_ML499M1-S1.tsv', sep = '\t', index_col=0)
# counts = pd.read_csv('../data/count_matrix_ML499M1-S1.tsv', sep = '\t', index_col=0)

outfile = ("../data/networks/ML499M1-S1_count_network.tsv")
cell_list_out = ("../data/networks/ML499M1-S1_cell_list.txt")

# number of highly variable genes to use before PCA
top_n = 1500

# number of latent dimensions for PCA
pca_dims = 50

# distance metric ("pearson", cosine)
dist_metric = "pearson"

# If proportion is not none, will prune using percentage (as in top k percent)
# Otherwise, uses all distances above given number
distance_use_proportion = 10
distance_use_cutoff = 0.5

################################################

# Output gene list
with open(cell_list_out, 'w') as f:
    for cell in cell_by_gene.index :
        f.write(f"{cell}\n")

# Clear duplicates
print("Pruning duplicates")
duplicate_mask = cell_by_gene.T.duplicated()
print("Number of duplicate columns:", duplicate_mask.sum())
cell_by_gene = cell_by_gene.loc[:, ~duplicate_mask]

# Get top genes
print("Getting highly variable genes")
gene_var = pd.Series(cell_by_gene.var(axis=0), index=cell_by_gene.columns)
top_genes = gene_var.sort_values(ascending=False).head(top_n).index
cell_by_gene_small = cell_by_gene[top_genes]

# Dimensionality reduction
print("Running dimensionality reduction")
pca = PCA(n_components=pca_dims)
cell_by_gene_pca = pca.fit_transform(cell_by_gene_small)
pc_names = [f"PC{i+1}" for i in range(cell_by_gene_pca.shape[1])]
cbg_df_pca = pd.DataFrame(cell_by_gene_pca, columns=pc_names, index=cell_by_gene_small.index)

print("Creating distance matrix")
if dist_metric == "pearson" :
    dist_matrix = cbg_df_pca.T.corr(method=dist_metric)

elif dist_metric == "cosine" :
    dist_array = cosine_similarity(cbg_df_pca)
    dist_matrix = pd.DataFrame(
        dist_array,
        index=cbg_df_pca.index,      # original cell names
        columns=cbg_df_pca.index
    )

if distance_use_proportion is not None :
    print("Calculating threshold")
    # exclude diagonal
    mask = ~np.eye(dist_matrix.shape[0], dtype=bool)
    off_diag_values = dist_matrix.values[mask]
    print(dist_matrix.shape)
    print(len(off_diag_values))

    k = int(len(off_diag_values) * (100 - distance_use_proportion) / 100)
    threshold = np.partition(off_diag_values, k)[k]
    if threshold < 0 :
        threshold = 0
    print(f"Threshold: {threshold}")

else :
    threshold = distance_use_cutoff

# Output
print("Outputting")
with open(outfile, "w") as out:
    for i in range(dist_matrix.shape[0]):
        for j in range(i + 1, dist_matrix.shape[0]):
            if dist_matrix.iloc[i, j] > threshold:
                cell_i = dist_matrix.index[i]
                cell_j = dist_matrix.columns[j]
                out.write(f"{cell_i}\t{cell_j}\t{dist_matrix.iloc[i, j]}\n")
