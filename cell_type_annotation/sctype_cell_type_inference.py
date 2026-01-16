#!/usr/bin/env python3

"""
Cell type identification using the sc-type algorithm
"""

### ---------------------------------------- ###

class sctype_like:
    
    """
    Class for automatic cell assignment given a group of positive and negative markers for each desired cell type
    Based on sc-type, https://www.nature.com/articles/s41467-022-28803-w
    
    Parameters
    ----------
    counts : scypy sparse csr_matrix
        Sparse matrix of cells x genes
    genes : list[str]
        List of gene identifiers as they appear in counts
    markers : pandas dataframe
        Pandas data frame with the following columns
        'tissue' : name of the tissue of origin
        'cell_type' : name of the cell type
        'positive_markers' : comma-separated string of positive markers
        'negative_markers' : comma-separated string of negative markers
        See https://sctype.app/database.php for sc-type curated database
    tissue : str
        Name of the tissue of interest
    """
    
    def __init__(self, counts, genes, markers, tissue):
        
        # Process markers
        
        tissue_markers, markers_specificity = self.parse_markers(markers, tissue, genes)
        
        self.tissue_markers = tissue_markers
        
        self.markers_specificity = markers_specificity
        
        # Process counts
        
        weighted_counts = self.parse_counts(counts, genes, markers_specificity)
        
        self.weighted_counts = weighted_counts
    
    ### ------------------------------------ ###
    ### PARSE MARKERS TABLE                  ###
    ### ------------------------------------ ###
    
    @staticmethod
    def parse_markers(mrks, t, detected_gns):
        
        # Replace nan values
        
        mrks = mrks.fillna('NONE')
        
        # Subset by tissue
        
        mrks_sub = mrks.loc[mrks['tissue'] == t,]
        
        # Get list of all markers (and implicit info about their occurrence)
        
        all_markers = (','.join(mrks_sub.loc[mrks_sub['positive_markers'] != 'NONE', 'positive_markers'].values) +
                       ',' +
                       ','.join(mrks_sub.loc[mrks_sub['negative_markers'] != 'NONE', 'negative_markers'].values))
        
        all_markers = [m for m in all_markers.split(',') if m in detected_gns]
        
        # Compute specificity
        
        occurrences = {m : all_markers.count(m) for m in set(all_markers)}
        
        min_o, max_o = min(occurrences.values()), max(occurrences.values())
        
        specificity = {m : 1 - (m_o - min_o) / (max_o - min_o) for m,m_o in occurrences.items()}
        
        return mrks_sub, specificity
    
    ### ------------------------------------ ###
    ### PARSE COUNT MATRIX                   ###
    ### ------------------------------------ ###
    
    def parse_counts(self, cnts, gns, m_specificity):
        
        # Subset counts matrix
        
        markers_idx = [gns.index(m) for m in m_specificity.keys()]
        
        cnts = cnts[:, markers_idx]
        
        # Scale counts matrix
        # N.B. will return a full matrix
        
        scaled_cnts = self.scale_features(cnts)
        
        # Annotate scaled_cnts
        
        scaled_cnts = pd.DataFrame(scaled_cnts, columns=m_specificity.keys())
        
        # Weight gene expression
        
        weighted_cnts = scaled_cnts * pd.Series(m_specificity)
        
        return weighted_cnts

    ### ------------------------------------ ###
    ### ASSIGN CELL TYPES                    ###
    ### ------------------------------------ ###
    
    def score_cell_identity(self):
        
        # Init score matrix
        
        scores = pd.DataFrame(np.zeros((self.weighted_counts.shape[0], self.tissue_markers.shape[0])), columns=self.tissue_markers['cell_type'].values)
        
        # Score individual cell types
        
        for _,(_,cell_type,pos,neg) in self.tissue_markers.iterrows():
            
            # Score positive markers
            
            pos = [p for p in pos.split(',') if p in self.markers_specificity.keys()]
            
            if not len(pos):
                
                pos_score = np.zeros(scores.shape[0])
            
            else:
                
                pos_score = self.weighted_counts[pos].sum(axis=1).values / np.sqrt(len(pos))
            
            # Score negative markers
            
            neg = [n for n in neg.split(',') if n in self.markers_specificity.keys()]
            
            if not len(neg):
                
                neg_score = np.zeros(scores.shape[0])
            
            else:
            
                neg_score = self.weighted_counts[neg].sum(axis=1).values / np.sqrt(len(neg))
            
            # Final score
            
            cell_type_score = pos_score - neg_score
            
            scores.loc[:, cell_type] = cell_type_score.copy()
    
        self.cell_scores = scores
    
        # Score identity based on highest value
    
        self.cell_identities = scores.idxmax(axis=1).values
    
    ### ------------------------------------ ###
    
    def reduce_to_cluster_level(self, clusters):
        
        # Merge cell identities scores and clusters
        
        cluster_scores = self.cell_scores.copy()
        
        cluster_scores.loc[:, 'cluster'] = clusters
        
        cluster_scores = cluster_scores.groupby('cluster').sum().reset_index(drop=False)
        
        self.cluster_scores = cluster_scores
        
        # Score identity based on highest value
        
        cluster_identities = pd.DataFrame({'cluster' : cluster_scores['cluster'].values,
                                           'identity' : cluster_scores.iloc[:, 1:].idxmax(axis=1).values})
        
        self.cluster_identities = cluster_identities
    
    ### ------------------------------------ ###
    ### UTILS                                ###
    ### ------------------------------------ ###
    
    def scale_features(self, mtx, gene_mean=[], gene_std=[]):
        
        # Z-score scaling of a sparse matrix
        # Note that after scaling the matrix will be dense since 0 values will be transformed to non-0
        
        # Get gene mean and std, if not provided
        if not len(gene_mean) or not len(gene_std):
            
            gene_mean, gene_std = self.get_gene_mean_and_std(mtx)
        
        # Scale features
        mtx_scaled = (mtx.toarray() - gene_mean) / (gene_std + 1e-6)
        
        return mtx_scaled
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_gene_mean_and_std(mtx):
        
        # Transposing matrix from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_matrix slicing is efficient along rows
        mtx = mtx.transpose()
        
        # Calculate features mean and variance
        # Remember that variance = mean of squared values - square of mean values
        mtx_squared = mtx.copy()
        mtx_squared.data **= 2
        gene_mean = np.asarray(mtx.mean(axis=1)).reshape(-1)
        gene_std = (np.asarray(mtx_squared.mean(axis=1)).reshape(-1) - (gene_mean ** 2))**0.5
        
        # Reshaping mean and std
        gene_mean = gene_mean.reshape((1, gene_mean.shape[0]))
        gene_std = gene_std.reshape((1, gene_std.shape[0]))
        
        return gene_mean, gene_std

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
