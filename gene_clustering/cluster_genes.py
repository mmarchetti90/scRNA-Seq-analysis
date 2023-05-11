#!/usr/bin/env python3

"""
Using dimensionality reduction and Leiden algorithm to find genes with similar expression,
similarly to how cells are clustered.

N.B. The matrix fed to the cluster_genes class must contain normalized counts.

N.B. The matrix fed to the cluster_genes class must be a genes x cells data frame, with each row a
gene and each column a cell. Starting columns reporting the gene id or name (can have both columns)
can be added. These starting columns must have the names 'GeneID' or 'GeneSymbol'
"""

### ---------------------------------------- ###

class cluster_genes:

    def __init__(self, gene_expression_matrix=None):
        
        try:
            
            counter = 0
        
            if 'GeneID' in gene_expression_matrix.columns:
                
                gene_ids = gene_expression_matrix.GeneID.values
                counter += 1
            
            else:
                
                gene_ids = [f'gene_{i}' for i in range(gene_expression_matrix.shape[0])]
                
            if 'GeneSymbol' in gene_expression_matrix.columns:
                
                gene_symbols = gene_expression_matrix.GeneSymbol.values
                counter += 1
            
            else:
                
                gene_symbols = [f'gene_{i}' for i in range(gene_expression_matrix.shape[0])]
            
            # Init data table
            self.analysis_data = pd.DataFrame({'GeneID' : gene_ids,
                                               'GeneSymbol' : gene_symbols})
            
            # Extract counts
            self.gene_expression_matrix = gene_expression_matrix.iloc[:, counter:].to_numpy()
        
        except:
            
            pass

    ### ------------------------------------ ###
    ### SAVE/LOAD RESULTS                    ###
    ### ------------------------------------ ###
    
    def load_data(self, file_dir='./'):
        
        # Load analysis data
        self.analysis_data = pd.read_csv(f'{file_dir}/analysis_data.tsv', sep='\t')
        
        # Load expression matrix
        self.gene_expression_matrix = pd.read_csv(f'{file_dir}/gene_expression_matrix.tsv.gz', header=None, index_col=None, sep='\t').to_numpy()
        
        # Load PCA info
        self.optimal_pca_components = int(open(f'{file_dir}/pca_info.tsv').read().split('\t')[-1])
        
        # Load PCA and UMAP models
        try:
            
            self.pca_model = pk.load(open(f'{file_dir}/pca.pkl', "rb"))
            self.umap_model = pk.load(open(f'{file_dir}/umap.pkl', "rb"))
        
        except:
            
            pass
        
        # Transform data
        self.reduce_dimensions(min_pca_components=20, max_pca_components=50, neighbors=30)
    
    ### ------------------------------------ ###
    
    def save_data(self, out_dir):
        
        if not exists(out_dir):
            
            mkdir(out_dir)
        
        # Save analysis data
        self.analysis_data.to_csv(f'{out_dir}/analysis_data.tsv', index=False, sep='\t')
        
        # Save expression matrix
        pd.DataFrame(self.gene_expression_matrix).to_csv(f'{out_dir}/gene_expression_matrix.tsv.gz', header=False, index=False, sep='\t')
        
        # Save PCA optimal components
        with open(f'{out_dir}/pca_info.tsv', 'w') as out:
            
            out.write(self.pca_info)
        
        # Save PCA and UMAP embeddings
        try:
            
            pk.dump(self.pca_model, open(f'{out_dir}/pca.pkl', "wb"))
            
        except:
            
            pass
        
        try:
            
            pk.dump(self.umap_model, open(f'{out_dir}/umap.pkl', "wb"))
                
        except:
            
            pass
    
    ### ------------------------------------ ###
    ### DIMENSIONALITY REDUCTION             ###
    ### ------------------------------------ ###
    
    def reduce_dimensions(self, min_pca_components=20, max_pca_components=50, neighbors=30):
        
        pca_input = self.gene_expression_matrix.copy()
        
        try:
            
            self.pca = self.pca_model.transform(pca_input)
            
        except:
            
            self.pca_model = PCA(max_pca_components)
            self.pca = self.pca_model.fit_transform(pca_input)
        
            # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
            x0, x1 = 0, min(len(self.pca_model.explained_variance_), max_pca_components - 1)
            y0, y1 = self.pca_model.explained_variance_[x0], self.pca_model.explained_variance_[x1]
            gradient = (y1 - y0) / (x1 - x0)
            intercept = y0
            difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(self.pca_model.explained_variance_[:max_pca_components])]
            self.optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
        
            # Elbow plot of explained variance
            plt.figure(figsize=(5, 3))
            plt.plot(range(x0 + 1, x1 + 2), self.pca_model.explained_variance_[:50] / 100, 'b', marker='o', markersize=5, linewidth=1)
            plt.plot([self.optimal_pca_components, self.optimal_pca_components], [y0 / 100, y1 / 100], linestyle='dashed', color='red', linewidth=1)
            plt.xlabel('PC')
            plt.ylabel('Explained Variance (%)')
            plt.tight_layout()
            plt.savefig("PCA_ExplainedVariance", dpi=300)
            plt.close()
            
            if self.optimal_pca_components < min_pca_components:
                
                self.optimal_pca_components = min_pca_components
        
        # Running UMAP with X neighbors
        try:
            
            self.umap = self.umap_models.transform(self.pca[:, :self.optimal_pca_components])
            self.umap = pd.DataFrame(data = self.umap, index = self.analysis_data.GeneSymbol.to_list(), columns = ['UMAP_1', 'UMAP_2'])
        
        except:
            
            self.umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
            self.umap = self.umap_model.fit_transform(self.pca[:, :self.optimal_pca_components])
            self.umap = pd.DataFrame(data = self.umap, index = self.analysis_data.GeneSymbol.to_list(), columns = ['UMAP_1', 'UMAP_2'])
        
        # Storing PCA info
        self.pca_info = f'optimal_pca_components\t{self.optimal_pca_components}'
    
    ### ------------------------------------ ###
    ### CLUSTERING                           ###
    ### ------------------------------------ ###
    
    def cluster_genes(self, n_neighbors=10):
        
        # Setting random seed (helps with clustering consistency)
        random.seed(42)
        
        # Computing kneighbors sparse matrix
        kneighbors_matrix = kneighbors_graph(X=self.pca[:, :self.optimal_pca_components], n_neighbors=n_neighbors)
        
        # Creating edges list
        sources, targets = kneighbors_matrix.nonzero()
        edges = list(zip(sources, targets))
        
        # Building igraph object
        graph = igraph.Graph(directed=True)
        graph.add_vertices(kneighbors_matrix.shape[0])
        graph.add_edges(edges)
        
        # Converting graph to undirected
        graph.to_undirected(mode='collapse', combine_edges=None)
        
        # Clustering using Leiden algorithm
        clusters = graph.community_leiden(objective_function='modularity', weights=None, resolution_parameter=1.0, beta=0.01, initial_membership=None, n_iterations=2, node_weights=None).membership
        
        # Updating analysis data
        self.analysis_data['Cluster'] = clusters
    
    ### ------------------------------------ ###
    ### PLOTTING                             ###
    ### ------------------------------------ ###

    def plot_clusters(self, save_dir='./', gene_ids=[], gene_symbols=[], clusters=[], dot_size=1.5):
        
        # Creating a filter for cells of interest
        if not len(gene_ids):
        
            gene_ids = self.analysis_data.GeneID
        
        if not len(gene_symbols):
        
            gene_symbols = self.analysis_data.GeneSymbol
        
        if not len(clusters):
        
            clusters = list(self.analysis_data.Cluster.unique())
            clusters.sort()
        
        cell_filter = list(self.analysis_data.GeneID.isin(gene_ids) & self.analysis_data.GeneSymbol.isin(gene_symbols) & self.analysis_data.Cluster.isin(clusters))

        # Splitting legend into legend_cols columns
        legend_cols = int(len(clusters) / 16) + 1

        # Adding expression data, then sorting by smallest value
        plot_data = self.umap.copy()
        plot_data['Clusters'] = list(self.analysis_data.Cluster)
    
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter,]
    
        # Finding cluster centers for writing text
        centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
    
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Clusters', palette='tab10', marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Clusters', ncol=legend_cols)
        for c_num,cl in enumerate(centers):
            x, y = cl
            plt.text(x, y, str(clusters[c_num]), horizontalalignment='center', size='small', color='black', weight='semibold')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{save_dir}/Clusters_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()

### ------------------MAIN------------------ ###

import igraph
import numpy as np
import pandas as pd
import pickle as pk
import random
import seaborn
import umap

from matplotlib import pyplot as plt
from os import mkdir
from os.path import exists
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph