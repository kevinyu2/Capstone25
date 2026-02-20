import scanpy as sc
import celltypist
import anndata as ad
import pandas as pd
import sys
import argparse
import os 
from sklearn.preprocessing import MinMaxScaler

parser = argparse.ArgumentParser(description="Process count and expression data.")

parser.add_argument(
    "--counts",
    type=str,
    required=True,
    help="Path to ATAC tsv"
)

parser.add_argument(
    "--expression",
    type=str,
    required=True,
    help="Path to expression tsv"
)

parser.add_argument(
    "--out-dir",
    type=str,
    required=True,
    help="Path to output directory"
)

parser.add_argument(
    "--model",
    type=str,
    default="Immune_All_Low.pkl",
    help="Name of cell type model"
)

parser.add_argument(
    "--num-top",
    type=int,
    default=1500,
    help="Name of (non-duplicate) genes to use"
)

args = parser.parse_args()

'''
Get cell annotations
'''
def get_cell_annotations(args, expression_df) :
    print("Getting cell annotations")
    # Download models if not done before
    model_dir = os.path.expanduser("~/.celltypist/models/")
    if not os.path.exists(model_dir) or not os.listdir(model_dir):
        print("No models found, downloading now")
        celltypist.models.download_models()

    # Turn to anndata
    adata = ad.AnnData(expression_df)
    adata.var_names = expression_df.columns
    adata.obs_names = expression_df.index

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Predict using model
    predictions = celltypist.annotate(
        adata,
        model=args.model,
        majority_voting=True
    )

    # Add predictions
    adata = predictions.to_adata()

    export_df = pd.DataFrame({
        'cell': adata.obs_names,
        'annotation': adata.obs['majority_voting']
    })

    # Save to CSV
    export_df.to_csv(f'{args.out_dir}/cell_annotations.csv', index=False)

    label_dict = dict(zip(export_df['cell'], export_df['annotation']))

    return label_dict

def make_umap(args, ad, name) :
    sc.settings.figdir = args.out_dir

    # Create umap
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.scale(ad)
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
    export_df.to_csv(f'{args.out_dir}/UMAP_coords_{name}.csv', index=False)

    # Save figures
    sc.pl.umap(
        ad, 
        color='cell_type', 
        title=f'UMAP {name}', 
        legend_loc='right', 
        frameon=False,
        show=False,
        save=f'UMAP_{name}.png'
    )

def run(args) :
    os.makedirs(args.out_dir, exist_ok=True)
    print("Reading in data")
    expression_df = pd.read_csv(args.expression, sep = '\t', index_col=0)
    counts_df = pd.read_csv(args.counts, sep = '\t', index_col=0)

    label_dict = get_cell_annotations(args, expression_df)

    # Choose only top expression genes
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(expression_df)
    gene_var = pd.Series(X_scaled.var(axis=0), index=expression_df.columns)
    top_n = args.num_top
    top_genes = gene_var.sort_values(ascending=False).head(top_n).index
    expression_var = expression_df[top_genes]

    # For top counts, first remove duplicates due to windows
    duplicate_mask = counts_df.T.duplicated()
    counts_unique = counts_df.loc[:, ~duplicate_mask]
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(counts_unique)
    gene_var = pd.Series(X_scaled.var(axis=0), index=counts_unique.columns)
    top_genes = gene_var.sort_values(ascending=False).head(top_n).index
    counts_var = counts_unique[top_genes]

    # Turn to anndatas
    expression_ad = ad.AnnData(expression_var)
    expression_ad.obs['cell_type'] = [label_dict[cell] for cell in expression_df.index]
    count_ad = ad.AnnData(counts_var)
    count_ad.obs['cell_type'] = [label_dict[cell] for cell in counts_df.index]

    make_umap(args, count_ad, "CopyNumber")
    make_umap(args, expression_ad, "Expression")
run(args)