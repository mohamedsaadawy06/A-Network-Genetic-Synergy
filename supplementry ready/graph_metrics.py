import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from data_loader import ROI_LABELS

def build_weighted_graph(matrix: np.ndarray, roi_labels: list) -> nx.Graph:
    """
    Build a weighted undirected graph from the connectivity matrix.
    Log-transform streamline counts to reduce skew.
    """
    W = np.log1p(matrix)
    G = nx.from_numpy_array(W)
    
    # Rename nodes to ROI labels
    mapping = {i: roi_labels[i] for i in range(len(roi_labels))}
    G = nx.relabel_nodes(G, mapping)
    return G

def compute_betweenness_centrality(G: nx.Graph) -> dict:
    """
    Compute Betweenness Centrality using inverse edge weights
    so high-streamline tracts act as short distances.
    """
    # Create distance attribute: distance = 1 / weight
    for u, v, d in G.edges(data=True):
        d['distance'] = 1.0 / d['weight'] if d['weight'] > 0 else 9999.0
        
    bc = nx.betweenness_centrality(G, weight='distance', normalized=True)
    return bc

def compute_effective_distance(matrix: np.ndarray,
                               roi_labels: list,
                               seed_label: str) -> dict:
    """
    Effective distance: D_eff(r, seed) = -log(P(seed -> r))
    where P is the row-normalized transition matrix.
    """
    # Row-normalize to transition probability matrix
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    P = matrix / row_sums

    # D_eff = -log(P_ij), treat zeros as infinity
    with np.errstate(divide='ignore'):
        D_eff_mat = -np.log(P)
    D_eff_mat[np.isinf(D_eff_mat)] = 9999.0

    # Shortest paths using Dijkstra
    seed_idx = roi_labels.index(seed_label)
    dist_from_seed = shortest_path(
        csr_matrix(D_eff_mat),
        method='D',
        indices=seed_idx,
        directed=False
    )
    return {roi_labels[i]: dist_from_seed[i] for i in range(len(roi_labels))}
