import scanpy as sc
import celltypist
import anndata as ad
import pandas as pd
import sys
import argparse
import os 
from sklearn.preprocessing import MinMaxScaler

# parser = argparse.ArgumentParser(description="Umap for repel.")

# parser.add_argument(
#     "--embeds",
#     type=str,
#     required=True,
#     help="Path to expression tsv"
# )

# parser.add_argument(
#     "--out-dir",
#     type=str,
#     required=True,
#     help="Path to output directory"
# )

# parser.add_argument(
#     "--model",
#     type=str,
#     default="Immune_All_Low.pkl",
#     help="Name of cell type model"
# )

# parser.add_argument(
#     "--num-top",
#     type=int,
#     default=1500,
#     help="Name of (non-duplicate) genes to use"
# )

# args = parser.parse_args()

sample = "HT268B1-Th1H3"
embeds = f"../data/embeds/{sample}_150_10_embed.tsv"
cell_annotations = f"./{sample}_UMAPs/cell_annotations.csv"
out_dir = f'./{sample}_UMAPs/'

'''
Get cell annotations
'''
def get_cell_annotations(cell_annotations) :
    label_dict = {}
    with open(cell_annotations, 'r') as f :
        for i,line in enumerate(f) :
            name, cell_type = line.rstrip().split(',')
            if i != 0 :
                label_dict[name] = cell_type

    return label_dict

def make_umap(ad, out_dir, name = "repel") :
    sc.settings.figdir = out_dir

    # Create umap
    sc.tl.pca(ad, svd_solver='arpack')
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=50)
    sc.tl.umap(ad)

    # Save coordinates
    ad.obs[['UMAP1', 'UMAP2']] = pd.DataFrame(ad.obsm['X_umap'], index=ad.obs_names)
    export_df = pd.DataFrame({
        'barcode': ad.obs_names,
        'cell_type': ad.obs['cell_type'],
        'UMAP1': ad.obs['UMAP1'],
        'UMAP2': ad.obs['UMAP2']
    })
    export_df.to_csv(f'{out_dir}/UMAP_coords_{name}.csv', index=False)

    # Save figures
    sc.pl.umap(
        ad,
        color='cell_type',
        title=f'UMAP {name}',
        legend_loc='right margin',
        frameon=False,
        show=False,
        save=f'UMAP_{name}.png'
    )

def run(embeds, cell_annotations,out_dir) :
    embed_df = pd.read_csv(embeds, sep = '\t', index_col=0)


    label_dict = get_cell_annotations(cell_annotations)

    # Turn to anndatas
    embed_ad = ad.AnnData(embed_df)
    embed_ad.obs['cell_type'] = [label_dict[cell] for cell in embed_df.index]
   
    make_umap(embed_ad, out_dir)
run(embeds, cell_annotations, out_dir)




    # # Choose only top expression genes
    # scaler = MinMaxScaler()
    # X_scaled = scaler.fit_transform(expression_df)
    # gene_var = pd.Series(X_scaled.var(axis=0), index=expression_df.columns)
    # top_n = args.num_top
    # top_genes = gene_var.sort_values(ascending=False).head(top_n).index
    # expression_var = expression_df[top_genes]

    # # For top counts, first remove duplicates due to windows
    # duplicate_mask = counts_df.T.duplicated()
    # counts_unique = counts_df.loc[:, ~duplicate_mask]
    # scaler = MinMaxScaler()
    # X_scaled = scaler.fit_transform(counts_unique)
    # gene_var = pd.Series(X_scaled.var(axis=0), index=counts_unique.columns)
    # top_genes = gene_var.sort_values(ascending=False).head(top_n).index
    # counts_var = counts_unique[top_genes]