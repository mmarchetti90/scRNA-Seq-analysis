#!/usr/bin/env python3

"""
Class for detection of a specific cell type at a cluster level
"""

### ---------------------------------------- ###

class find_cell_type:
    
    def __init__(self, expression_matrix, gene_names, clusters, gene_signature=[], reference_expression=[]):
        
        self.expression_matrix = expression_matrix
        self.gene_names = gene_names
        self.clusters = clusters
        
        self.gene_signature = gene_signature
        
        if type(reference_expression) != pd.core.frame.DataFrame:
            
            reference_expression = pd.DataFrame(reference_expression)
            
        self.reference_expression = reference_expression
        
        if reference_expression.shape[0] > 0:
            
            # Mapping genes in expression_matrix to genes in reference_expression
            self.gene_map = [gene_names.index(gene) for gene in reference_expression.GeneName.values]
    
    ### ------------------------------------ ###
    
    def analyze_data(self, gene_pool=[], bins=50, ctrl_genes_num=50, max_cells_in_memory=100000, top_fraction=0.33, correl_p_thr=0.05):
        
        # Init results dictionary
        self.results = {'cluster_id' : np.sort(np.unique(self.clusters))}
        
        ### gene_signature dependent methods
        
        if len(self.gene_signature) > 0:
            
            # Score cells for gene_signature
            print('Scoring gene signature')
            signature_score = self.score_gene_set(gene_pool, bins, ctrl_genes_num, max_cells_in_memory)
            self.results['gene_signature'] = signature_score
            
            # AUC of gene_signature recovery curve
            print('Scoring AUC')
            auc_scores = self.auc_gene_set(top_fraction = 0.33)
            self.results['auc_scores'] = auc_scores
        
        ### reference_expression dependent metrics
        
        if self.reference_expression.shape[0] > 0:
            
            # Score cells via euclidean distance to reference_expression
            print('Scoring distance from reference')
            distance_from_ref = self.euclidean_distance_from_reference()
            self.results['ref_distance'] = distance_from_ref
            
            # Score cells via Pearson correlation to reference_expression
            print('Scoring Pearson correlation with reference')
            pearson_r, pearson_r2, pearson_pval, pearson_padj = self.pearson_correlation_to_reference()
            self.results['pearson_r'] = pearson_r
            self.results['pearson_r2'] = pearson_r2
            self.results['pearson_pval'] = pearson_pval
            self.results['pearson_padj'] = pearson_padj
            
            # Score cells via Spearman correlation to reference_expression
            print('Scoring Spearman correlation with reference')
            spearman_r, spearman_r2, spearman_pval, spearman_padj = self.spearman_correlation_to_reference()
            self.results['spearman_r'] = spearman_r
            self.results['spearman_r2'] = spearman_r2
            self.results['spearman_pval'] = spearman_pval
            self.results['spearman_padj'] = spearman_padj
        
        self.results = pd.DataFrame(self.results)
        
        # Compound scores
        self.compund_scores(correl_p_thr)
        
        return self.results
    
    ### ------------------------------------ ###
    ### GENE SET SCORE                       ###
    ### ------------------------------------ ###
    
    def score_gene_set(self, gene_pool=[], bins=50, ctrl_genes_num=50, max_cells_in_memory=100000):

        # Subsetting gene_set and gene_list for detected genes
        gene_set = [gene for gene in self.gene_signature if gene in self.gene_names]
        gene_pool = [gene for gene in gene_pool if gene in self.gene_names]

        # Making sure that gene_set and gene_pool have len > 0
        if not len(gene_set):
            
            print("Error: empty gene_set")
            return None
        
        if not len(gene_pool):
            
            gene_pool = self.gene_names
        
        # Mean of each scaled gene in the gene_pool across all cells
        # Slow implementation, but manageable when theres a lot of cells
        if self.expression_matrix.shape[0] <= max_cells_in_memory:
            
            gene_means = self.scale_features(self.expression_matrix).mean(axis=0)
        
        else:
            
            # N.B. Could have used the .mean() method, but for some reason it kept returning the wrong results (so did the .sum() method)
            #gene_means = np.array([self.scale_features(self.expression_matrix[:, self.gene_names.index(g)]).mean(axis=0)[0] for g in gene_pool])
            gene_means = np.ravel(np.array([sum(self.scale_features(self.expression_matrix[:, self.gene_names.index(g)])) / self.expression_matrix.shape[0] for g in gene_pool]))
        
        # Rank genes based on binned expression level, then for each bin of the genes in the gene_set, pick ctrl_genes random genes for genes with matched binned expression
        bin_size = len(gene_means) // bins
        expression_order = np.argsort(gene_means)
        
        set_order_indexes = []
        for gene in gene_set:
            
            index = expression_order[self.gene_names.index(gene)]
            
            set_order_indexes.append(index)
        
        ctrl_indexes = []
        for index in set_order_indexes:
            
            # Fix this. It's a numpy array, not a pandas df
            random_pool = np.where((expression_order >= (index - bin_size // 2)) &
                                   (expression_order <= (index + bin_size // 2)) &
                                   ~(np.isin(expression_order, set_order_indexes)))[0].tolist()
            
            np.random.shuffle(random_pool)
            ctrl_indexes.extend(random_pool[:ctrl_genes_num])
        ctrl_indexes = list(set(ctrl_indexes)) # Removing duplicates
        
        # Computing the mean of gene_set genes for each cell
        set_indexes = [self.gene_names.index(g) for g in gene_set if g in self.gene_names]
        if self.expression_matrix.shape[0] <= max_cells_in_memory:
            
            set_means = self.scale_features(self.expression_matrix[:, set_indexes]).mean(axis=1)

        else:

            set_means = np.array([self.scale_features(self.expression_matrix[:, g]).reshape(-1) for g in set_indexes]).T.mean(axis=1)

        # Computing the mean of ctrl_set genes for each cell
        if self.expression_matrix.shape[0] <= max_cells_in_memory:
            
            ctrl_means = self.scale_features(self.expression_matrix[:, ctrl_indexes]).mean(axis=1)

        else:

            ctrl_means = np.array([self.scale_features(self.expression_matrix[:, g]).reshape(-1) for g in ctrl_indexes]).T.mean(axis=1)
        
        set_score = set_means - ctrl_means
        
        # Computing cluster level values
        cl_set_score = []
        for cl in np.sort(np.unique(self.clusters)):
            
            cl_set_score.append(np.median(set_score[self.clusters == cl]))
        
        return np.array(cl_set_score)
    
    ### ------------------------------------ ###
    
    def scale_features(self, mtx, gene_mean=[], gene_std=[]):
        
        # Z-score scaling of a sparse matrix
        # Note that after scaling the matrix will be dense since 0 values will be transformed to non-0
        
        # Get gene mean and std, if not provided
        if not len(gene_mean) or not len(gene_std):
            
            gene_mean, gene_std = self.get_gene_mean_and_std(mtx)
        
        # Scale features, parsing chunks of max_cells_in_memory
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
    
    ### ------------------------------------ ###
    ### GENE SET AUC                         ###
    ### ------------------------------------ ###
    
    def auc_gene_set(self, top_fraction = 0.33):
        
        max_genes = int(self.expression_matrix.shape[1] * top_fraction)
        auc_scores = []
        for cl in np.sort(np.unique(self.clusters)):
            
            # Extract cluster gene expression and rank them
            cl_counts = np.median(self.expression_matrix[self.clusters == cl, ].toarray(), axis=0)
            gene_rank = np.argsort(cl_counts)[::-1]
            
            # Randomize genes with same expression
            start = 0
            while start < max_genes - 1:
                
                stop = start + 1
                while cl_counts[gene_rank[start]] == cl_counts[gene_rank[stop]]:
                    
                    if stop == len(gene_rank) - 1:
                        
                        stop += 1
                        break
                    
                    stop += 1
            
                if stop != start + 1:
                    
                    gene_rank[start : stop] = np.random.choice(gene_rank[start : stop], size=(stop - start), replace=False)
                
                start = stop

            # Get top genes
            top_genes = np.array(self.gene_names)[np.argsort(cl_counts)[::-1][:max_genes]]
            
            # Calculate recovery curve of gene from the set
            recovery_curve = [sum([1 for gene in self.gene_signature if gene in top_genes[:i]]) for i in range(len(top_genes))]
            
            # Calculate Area Uder the Curve (AUC)
            cl_auc = sum(recovery_curve)

            auc_scores.append(cl_auc)

        # Normalize AUC to max value
        auc_scores = np.array(auc_scores) / max(auc_scores)
        
        return auc_scores
    
    ### ------------------------------------ ###
    ### EUCLIDEAN DISTANCE                   ###
    ### ------------------------------------ ###
    
    def euclidean_distance_from_reference(self):
        
        distances = []
        for cl in np.sort(np.unique(self.clusters)):
            
            cl_counts = np.median(self.expression_matrix[self.clusters == cl, ][:, self.gene_map].toarray(), axis=0)
        
            cl_distance = np.sqrt(np.sum(np.power(cl_counts - self.reference_expression.iloc[:, -1].values, 2)))
            
            distances.append(cl_distance)
        
        return np.array(distances)
    
    ### ------------------------------------ ###
    ### CORRELATION                          ###
    ### ------------------------------------ ###
    
    def pearson_correlation_to_reference(self):
        
        r_vals, r2_vals, p_vals = [], [], []
        for cl in np.sort(np.unique(self.clusters)):
            
            cl_counts = np.median(self.expression_matrix[self.clusters == cl, ][:, self.gene_map].toarray(), axis=0)
        
            pearson_r, pearson_pval = pearsonr(cl_counts, self.reference_expression.iloc[:, -1].values, alternative='two-sided')
            pearson_r2 = pearson_r**2
            
            r_vals.append(pearson_r)
            r2_vals.append(pearson_r2)
            p_vals.append(pearson_pval)
        
        _, p_adjs = fdrcorrection(p_vals, alpha=0.05, is_sorted=False)
        
        return np.array(r_vals), np.array(r2_vals), np.array(p_vals), np.array(p_adjs)
    
    ### ------------------------------------ ###
    
    def spearman_correlation_to_reference(self):
        
        r_vals, r2_vals, p_vals = [], [], []
        for cl in np.sort(np.unique(self.clusters)):
            
            cl_counts = np.median(self.expression_matrix[self.clusters == cl, ][:, self.gene_map].toarray(), axis=0)
        
            spearman_r, spearman_pval = spearmanr(cl_counts, self.reference_expression.iloc[:, -1].values, alternative='two-sided')
            spearman_r2 = spearman_r**2
            
            r_vals.append(spearman_r)
            r2_vals.append(spearman_r2)
            p_vals.append(spearman_pval)
        
        _, p_adjs = fdrcorrection(p_vals, alpha=0.05, is_sorted=False)
        
        return np.array(r_vals), np.array(r2_vals), np.array(p_vals), np.array(p_adjs)
    
    ### ------------------------------------ ###
    ### COMPOUND METRICS                     ###
    ### ------------------------------------ ###
    
    def compund_scores(self, correl_p_thr=0.05):
        
        # Combine scores after transforming with softmax/min
        compound_scores = np.ones(len(set(self.clusters)))
        for metric in ['gene_signature', 'auc_scores', 'ref_distance', 'pearson_r', 'spearman_r']:
        
            if metric not in self.results.columns:
                
                continue
            
            cluster_values = self.results.loc[:, metric].values

            if metric in ['pearson_r', 'spearman_r']:
                
                cluster_values[self.results.loc[:, metric.replace('_r', '_padj')].values > correl_p_thr] = 0.
            
            # Computing softmax (or softmin for ref_distance)
            if metric == 'ref_distance':
                
                transformed_values = [np.exp(- d) for d in cluster_values]
                scores = [tv / sum(transformed_values) for tv in transformed_values]
            
            else:
                
                transformed_values = [np.exp(d) for d in cluster_values]
                scores = [tv / sum(transformed_values) for tv in transformed_values]
            
            compound_scores *= scores
            
        final_rank = rankdata(- compound_scores)
        
        self.results['compound_scores'] = compound_scores
        self.results['final_rank'] = final_rank
    
    ### ------------------------------------ ###
    ### DIAGNOSTIC PLOTS                     ###
    ### ------------------------------------ ###
    
    def diagnostic_plots(self, correl_p_thr=0.05):
        
        # Define figsize
        figsize = (0.33 * self.results.shape[0], 0.1 * self.results.shape[0])
        
        # Define y labels
        y_labels = {'gene_signature' : 'Gene set signature',
                    'auc_scores' : 'Normalized AUC',
                    'ref_distance' : 'Distance from reference',
                    'pearson_r' : 'Pearson correlation to reference',
                    'spearman_r' : 'Spearman correlation to reference',
                    'compound_scores' : 'Compound cluster scores'}
        
        # Plotting
        for metric in y_labels.keys():
            
            if metric not in self.results.columns:
                
                continue
            
            y_vals = self.results.loc[:, metric].values
            y_lab = y_labels[metric]
            
            if metric in ['pearson_r', 'spearman_r']:
                
                y_vals[self.results.loc[:, metric.replace('_r', '_padj')].values > correl_p_thr] = 0.
            
            # Plot of compound scores
            plt.figure(figsize=figsize)
            plt.plot(self.results.cluster_id.values,
                     y_vals,
                     marker='o',
                     c='red',
                     markersize=6,
                     markeredgecolor='black',
                     markeredgewidth=1,
                     linewidth=0)
            plt.xlabel('Cluster')
            plt.ylabel(y_lab)
            plt.xticks(ticks=self.results.cluster_id.values.astype(int))
            plt.savefig(f'{metric}.png', bbox_inches='tight', dpi=300)
            plt.close()
            
            if metric == 'compound_scores':
                
                cluster_order = np.argsort(- self.results.compound_scores.values)

                # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
                x0, x1 = int(self.results.final_rank.min()) - 1, int(self.results.final_rank.max()) - 1
                y0, y1 = self.results.compound_scores.max() , self.results.compound_scores.min()
                gradient = (y1 - y0) / (x1 - x0)
                intercept = y0
                difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(self.results.loc[:, 'compound_scores'].values[cluster_order])]
                cluster_cutoff = difference_vector.index(max(difference_vector))
                print(f'Good clusters: {", ".join(cluster_order[:cluster_cutoff].astype(str))}')

                # Elbow plot of explained variance
                plt.figure(figsize=figsize)
                plt.plot(range(x0, x1 + 1), self.results.loc[:, 'compound_scores'].values[cluster_order], 'b', marker='o', markersize=5, linewidth=1)
                plt.plot([cluster_cutoff - 0.5, cluster_cutoff - 0.5], [y0, y1], linestyle='dashed', color='red', linewidth=1)
                plt.xlabel('Cluster')
                plt.ylabel('Compound cluster scores')
                plt.xticks(ticks=range(x0, x1 + 1), labels=cluster_order)
                plt.tight_layout()
                plt.savefig('compound_score_elbow_plot.png', dpi=300)
                plt.close()
        
### ------------------MAIN------------------ ###

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from scipy.stats import pearsonr, rankdata, spearmanr
from statsmodels.stats.multitest import fdrcorrection
