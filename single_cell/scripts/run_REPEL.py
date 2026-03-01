import numpy as np
import pandas as pd
import networkx as nx
import scipy
from scipy.spatial.distance import pdist, squareform
import sklearn
import matplotlib.pyplot as plt



def load_network_sparse(cell_index,net_file,ngene):
    ppi_df = pd.read_csv(net_file,header=None,sep='\t')
    ppi_df.columns = ["cell1", "cell2", "weight"]

    
    A = np.zeros((ngene,ngene))
    for i, row in ppi_df.iterrows() :
        row_idx = cell_index[row['cell1']] 
        col_idx = cell_index[row['cell2']]
        A[row_idx, col_idx] = row['weight']
        A[col_idx, row_idx] = row['weight']
    assert (A == A.T).all()
    zero_rows = np.all(A == 0, axis=1)
    diag_indices = np.arange(ngene)
    A[diag_indices[zero_rows], diag_indices[zero_rows]] = 1
    return A


def load_all_nets(cell_index,ppi_files,n_gene):
    '''
    parameters:
    - ppi_files: [str, str, ...], list of network file paths, each file should contain three columns: [protein1, protein2, score]
    - ref_gene_file: str, file path, the file contains all genes, one gene per line
    output:
    - nets: n_file x n_gene x n_gene array with ppi networks
    '''
    n_file = len(ppi_files)
    nets = np.zeros((n_file,n_gene,n_gene))
    for i in range(n_file):
        A = load_network_sparse(cell_index,ppi_files[i],n_gene)
        nets[i,:,:] = A
    return nets


def compute_rwr_original_sparse(ppi_files,restart_prob,ngene,nets):
    ''' 
    - ppi_files: list of network file paths
    - restart_prob: RWR restart probability
    - ngene: number of genes
    output:
    - walks: for the i-th RWR result, [i,:,:], each column is the stationary distribution of a node
    '''
    n_file = len(ppi_files)
    e = np.ones(ngene)
    I = np.eye(ngene)
    walks = np.zeros((n_file,ngene,ngene))
    for i in range(n_file):
        A = nets[i,:,:]
        d = A @ e
        P = A / d # transition matrix
        W = (I - (1 - restart_prob) * P)
        W = np.linalg.inv(W)
        W = W * restart_prob 
        walks[i,:,:] = W
    return walks

def svd_embed_sparse_func(walks, ngene, embed_dim):
    n_net = walks.shape[0]
    mat = np.zeros((ngene,ngene))
    W_updated = np.zeros_like(walks)
    for i in range(n_net):
        W = walks[i,:,:]
        W[W<=1e-8] = 0
        W = np.log(W, where = W > 1e-8)
        W_updated[i,:,:] = W
        tmp = W.T @ W
        mat = mat + tmp
    eigenvalues, eigenvectors = scipy.sparse.linalg.eigs(mat,k=embed_dim)
    x = np.diag(np.sqrt(np.sqrt(eigenvalues))) @ eigenvectors.T
    return np.real(x)



def augment_graph(nets, ngene, gene_clusters, mustlink_weight, cannotlink_weight):
    '''
    - nets: original adjacency matrices directly read from PPI files
    - gene_clusters: (num_clusers, num_genes), binary matrix indicating which gene belongs to which clusters
    '''
    n_nets = nets.shape[0]
    n_clusters = gene_clusters.shape[0]
    augmented = np.zeros((n_nets,(ngene+n_clusters),(ngene+n_clusters)))
    for i in range(n_nets):
        A = nets[i,:,:]
        A_block = np.block([[A,mustlink_weight*gene_clusters.T],[mustlink_weight*gene_clusters,cannotlink_weight*np.ones((n_clusters,n_clusters))]])
        np.fill_diagonal(A_block,0)
        zero_rows = np.all(np.absolute(A_block) == 0, axis=1)
        diag_indices = np.arange(ngene+n_clusters)
        A_block[diag_indices[zero_rows], diag_indices[zero_rows]] = 1
        augmented[i,:,:] = A_block
    return augmented

def augmented_RWR(augmented_nets, restart_prob):
    '''
    RWR for augmented graph which contains negative edge weights
    '''
    n_nets = augmented_nets.shape[0]
    n_nodes = augmented_nets.shape[1]
    augmented_walks = np.zeros((n_nets,n_nodes,n_nodes))
    e = np.ones(n_nodes)
    for i in range(n_nets):
        A = augmented_nets[i,:,:]
        d = np.absolute(A) @ e
        L = np.diag(d) - (1-restart_prob)*A
        L_inv = np.linalg.inv(L)
        W = restart_prob*(np.diag(d) @ L_inv)
        augmented_walks[i,:,:] = W
    return augmented_walks


def augmented_SVD_with_cannolink(aug_walks, embed_dim):
    n_net = aug_walks.shape[0]
    n_node = aug_walks.shape[1]
    mat = np.zeros((n_node,n_node))
    W_updated = np.zeros_like(aug_walks)
    for i in range(n_net):
        W = aug_walks[i,:,:]
        min_entry = W.min()
        if min_entry > 0:
            min_entry = 0.0
        W = W - min_entry
        W[W<=1e-8] = 0
        W = np.log(W, where = W > 1e-8)
        W_updated[i,:,:] = W
        tmp = W.T @ W
        mat = mat + tmp
    eigenvalues, eigenvectors = scipy.sparse.linalg.eigs(mat,k=embed_dim)
    x = np.diag(np.sqrt(np.sqrt(eigenvalues))) @ eigenvectors.T
    return np.real(x)



'''
for each augmented node, randomly choose a fixed number of nodes to connect to
the number of genes that each augmented node connects to are the same except for the last one
'''
def random_split_vector(n_gene, num_sub_vectors,seed=None):
    if seed is not None:
        np.random.seed(seed)
    
    input_vector = [i for i in range(n_gene)]
    
    if num_sub_vectors <= 0 or num_sub_vectors > len(input_vector):
        raise ValueError("Invalid number of sub-vectors")
    
    shuffled_vector = np.random.permutation(input_vector)
    sub_vector_size = len(shuffled_vector) // num_sub_vectors
    
    group_matrix = np.zeros((num_sub_vectors, len(input_vector)), dtype=int)
    res_matrix = np.zeros((num_sub_vectors, n_gene), dtype=int)
    
    start_index = 0
    for i in range(num_sub_vectors):
        end_index = start_index + sub_vector_size
        
        if i == num_sub_vectors - 1:
            end_index = len(shuffled_vector)
        
        selected_indices = shuffled_vector[start_index:end_index] 
        group_matrix[i, np.isin(shuffled_vector, selected_indices)] = 1
        
        start_index = end_index
    
    res_matrix[:,shuffled_vector] = group_matrix
    
    return res_matrix



# Just gets the embeddings per fold, returns as a list of embeddings for each fold
def run_dim_reduction_pipeline(cell_index,ppi_files,n_gene,method=None,restart_prob=None,embed_dim=None,n_cluster=None):
    ''' 
    parameters:
    - ppi_files: list of str, list of file paths to ppi networks
    - n_gene: int, number of genes
    - method: list of str, one or more of Mashup, REPEL
    - restart_prob: float, RWR restart probability
    - embed_dim: int, number of dimension
    - rand: int, random split
    - org: str, "yeast" or "Ecoli" 
    - n_fold: int, total number of folds
    - ont_type: str, bp or mf or cc
    - ont_size1: int, 11, 31, 101
    - ont_size2: int, 30, 100, 300
    - n_cluster: int, number of random augmented nodes
    output:
    - performance_dict: a dictionary contains list of performances for all methods
    '''


    if method == "Mashup":
        print("Mashup")
        nets = load_all_nets(cell_index,ppi_files,n_gene)
        walks = compute_rwr_original_sparse(ppi_files,restart_prob,n_gene,nets)
        x = svd_embed_sparse_func(walks, n_gene, embed_dim)
        embed = x
    elif method == "REPEL":
        print("REPEL")
        nets = load_all_nets(cell_index,ppi_files,n_gene)
        rand_cluster = random_split_vector(n_gene, n_cluster,seed=None)
        rand_graph = augment_graph(nets, n_gene, rand_cluster, 1, -1)
        rand_rwr_res = augmented_RWR(rand_graph, restart_prob)
        x = augmented_SVD_with_cannolink(rand_rwr_res, embed_dim)
        embed = x

    else:
        print("Haven't implemented yet")
        return
    

    return embed




def get_embeddings_repel(cell_index,ppi_files, n_gene, restart_prob=0.5, embed_dim=400, n_cluster=15):
    method= 'REPEL'

    embeds = run_dim_reduction_pipeline(cell_index,ppi_files,n_gene,method=method,restart_prob=restart_prob,embed_dim=embed_dim,n_cluster=n_cluster)

    embeds_no_augment = embeds[:, :n_gene]
    return embeds_no_augment




sample = "ML499M1-S1"
network_files = [f'../data/networks/{sample}_count_network.tsv', f'../data/networks/{sample}_expression_network.tsv']
cell_list_file = f'../data/networks/{sample}_cell_list.txt'
embed_dim = 150
n_cluster = 5
out_file = f'../data/embeds/{sample}_{embed_dim}_{n_cluster}_embed.tsv'


cell_index = {}
cell_list = []
with open(cell_list_file, 'r') as f :
    for i, line in enumerate(f) :
        cell_index[line.rstrip()] = i
        cell_list.append(line.rstrip())

n_gene = len(cell_index)


obj = get_embeddings_repel(cell_index,network_files, n_gene)
obj = obj.T
latent_cols = [f"latent_dim_{i}" for i in range(obj.shape[1])]
result_df = pd.DataFrame(obj, index=cell_list, columns=latent_cols)
result_df.to_csv(out_file, sep = '\t')
