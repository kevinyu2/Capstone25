Create matrices:

1. make_matrix.R (get cell by gene expression matrix from RDS)
2. cell_overlap.ipynb (choose only the cells in both count and expression data)

Create networks: 

3. network_create.py (use PCA and distances to create networks)
(optional) make_networks.ipynb (for exploration)
(optional) yeast_network_metrics.ipynb (stats for the original networks if we want to get similar patterns)

Run repel:

4. run_REPEL.py (create repel embedding)

Create UMAPs:

5. non_repel.py (get cell annotations from expression data and make umaps from just expression and count)
    Ex: python non_repel.py --counts ../data/count_matrix_HT268B1-Th1H3.tsv --expression ../data/expression_raw_matrix_HT268B1-Th1H3.tsv --out-dir ./HT268B1-Th1H3_UMAPs
6. repel_umap.py (make umap from repel embeddings)
