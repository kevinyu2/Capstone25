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
6. repel_umap.py (make umap from repel embeddings)
