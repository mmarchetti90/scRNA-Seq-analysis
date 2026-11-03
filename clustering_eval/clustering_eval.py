#!/usr/bin/env python3

"""
Evaluate n_neighbors and resolution parameters for Leiden clustering of single-cell RNA data
N.B. Uses metrics that don't take biology into consideration, so proposed metrics have to be evaluated
"""

### IMPORTS -------------------------------- ###

import igraph
import numpy as np
import pandas as pd
import random

from datetime import datetime
from itertools import combinations
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, silhouette_score
from sklearn.neighbors import kneighbors_graph

### CLASSES -------------------------------- ###

class LeidenClusteringEvaluator:

    def __init__(self) -> None:
        pass

    ### ------------------------------------ ###

    def evaluate_parameters(
        self,
        pca_data: pd.DataFrame,
        n_neighbors_vals: list[int]=[10, 20, 30, 40, 50],
        resolution_vals: list[float]=[0.2, 0.4, 0.6, 0.8, 1.0],
        subset_fraction: float=0.75,
        n_iterations: int=10
    ) -> pd.DataFrame:

        # Evaluate clustering
        self.log_trace('Starting clustering evaluation')
        metrics, metrics_header = [], ['n_neighbors', 'resolution', 'silhouette', 'ari', 'ami']
        subset_size = round(subset_fraction * pca_data.shape[0])
        for n_neighbors in n_neighbors_vals:
            for resolution in resolution_vals:
                self.log_trace(f'Evaluating for n_neighbors = {n_neighbors} & resolution = {resolution}')
                # Cluster using random subsets of cells
                clustering_results = []
                for _ in range(n_iterations):
                    indexes_sub = np.random.choice(pca_data.index.values, subset_size, False)
                    pca_data_sub = pca_data.loc[indexes_sub].copy()
                    new_clusters = self.cluster_cells(pca_data_sub, n_neighbors, resolution)
                    new_clusters = pd.Series(new_clusters, index=indexes_sub)
                    if len(set(new_clusters)) >= 2:
                        clustering_results.append(new_clusters)
                if not len(clustering_results):
                    continue
                clustering_results = pd.concat(clustering_results, axis=1).fillna(-1).astype(int)
                # Evaluate Silhouette
                silhouette = []
                for col in clustering_results.columns:
                    pca_data_sub = pca_data.loc[clustering_results.index[clustering_results[col] != -1]].copy()
                    clusters = clustering_results.loc[clustering_results[col] != -1, col].values
                    silhouette.append(silhouette_score(pca_data_sub, clusters, metric='euclidean', sample_size=10000))
                silhouette = np.median(silhouette)
                # Evaluate stability using following metrics
                # Median of pairwise Adjusted Rand Index (N.B. doesn't care if cluster count or labels match, just overlap of labels)
                # Median of pairwise Adjusted Mutual Info score
                ari, ami = [], []
                for a,b in combinations(clustering_results.columns, 2):
                    clustering_sub = clustering_results.loc[(clustering_results[a] != -1) & (clustering_results[b] != -1), [a, b]].copy()
                    ari.append(adjusted_rand_score(clustering_sub[a].values, clustering_sub[b].values))
                    ami.append(adjusted_mutual_info_score(clustering_sub[a].values, clustering_sub[b].values))
                ari, ami = np.median(ari), np.median(ami)
                # Store data
                metrics.append((n_neighbors, resolution, silhouette, ari, ami))
                self.log_trace(trace_message='  \u2714 Done')
        metrics = pd.DataFrame(metrics, columns=metrics_header)
        
        return metrics

    ### ------------------------------------ ###

    @staticmethod
    def cluster_cells(
        pca_data: pd.DataFrame,
        n_neighbors=10,
        resolution_parameter=1.0,
        beta=0.01
    ) -> list[int]:
        
        # Setting random seed (helps with clustering consistency)
        random.seed(42)
        
        # Computing kneighbors sparse matrix
        kneighbors_matrix = kneighbors_graph(X=pca_data, n_neighbors=n_neighbors)
        
        # Creating edges list
        sources, targets = kneighbors_matrix.nonzero() # type: ignore
        edges = list(zip(sources, targets))
        
        # Building igraph object
        graph = igraph.Graph(directed=True)
        graph.add_vertices(kneighbors_matrix.shape[0])
        graph.add_edges(edges)
        
        # Converting graph to undirected
        graph.to_undirected(mode='collapse', combine_edges=None)
        
        # Clustering using Leiden algorithm
        clusters = graph.community_leiden(
            objective_function='modularity',
            weights=None,
            resolution_parameter=resolution_parameter,
            beta=beta,
            initial_membership=None,
            n_iterations=2,
            node_weights=None
        ).membership

        return clusters

    ### ------------------------------------ ###
    
    @staticmethod
    def log_trace(trace_message: str, timestamp: bool=True) -> None:

        """
        Logs a time-stamped trace message to the console

        Parameters
        ----------
        trace_message: str
            The message to log
        timestamp: bool=True
            Set to true to prepend a timestamp
        """

        if timestamp:
            trace_timestamp = datetime.now().strftime("%H:%M:%S")
            timestamped_trace_message = f"[{trace_timestamp}] {trace_message}"
            print(timestamped_trace_message)
        else:
            print(trace_message)

### MAIN ----------------------------------- ###

if __name__ == '__main__':
    
    # Generate mock data
    pca_data = pd.DataFrame({f'PCA{i+1}' : np.random.normal(0., 1., 10000) for i in range(20)})
    pca_data.loc[np.random.choice([True, False], pca_data.shape[0]), :] += 1
    
    # Evaluate clustering
    evaluator = LeidenClusteringEvaluator()
    evaluation = evaluator.evaluate_parameters(pca_data)
    evaluation.to_csv('test_clustering_eval.tsv', sep='\t', index=False)
