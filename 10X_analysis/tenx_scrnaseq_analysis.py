#!/usr/bin/env python3

"""
DESCRIPTION:
Class to read in and analyze filtered HDF5 count file from Cellranger

NOTES:
HDF5 files from Cellranger have the following structure:
    - Dictionary with one key, 'matrix';
        - Dictionary with the following keys: 'barcodes', 'data', 'features', 'indices', 'indptr', 'shape';
            - 'barcodes' dataset of cell barcodes;
            - 'data' dataset of the counts for each feature of each cell (only non-zero values are reported);
            - 'features' dictionary with the following keys: '_all_tag_keys', 'feature_type', 'genome', 'id', 'name';
                - 'feature_type' dataset of the type of each feature (e.g. 'Gene Expression');
                - 'genome' dataset of the name of the genome reference for each feature;
                - 'id' dataset of feature IDs (e.g. Ensembl gene name);
                - 'name' dataset of feature names (e.g. gene name);
            - 'indices' dataset of the feature index related to each count reported in the 'data' dataset;
            - 'indptr' dataset of the index of the last 'data' associated with a cell (e.g. cell_1 counts = 'data'[0 : 3885]);
            - 'shape' dataset indicating the shape of the count matrix (i.e. features x barcodes);
"""

### ---------------------------------------- ###

class tenx_scranseq:
    
    def __init__(self):
        
        # Creating color palette for up to 60 clusters
        cluster_colors = [(r/3, g/3, b/3) for r in range(4) for g in range(4) for b in range(4)][1:-1]
        cluster_colors = [cluster_colors[i1 + i2] for i2 in range(4) for i1 in range(0, len(cluster_colors), 4) if i1 + i2 < len(cluster_colors)] # Spacing colors
        random_indexes = [33, 18, 41, 43, 34, 8, 56, 31, 16, 25, 2, 7, 46, 36, 5, 22, 11, 26, 53, 9, 14, 19, 58, 23, 1, 42, 55, 61, 35, 49, 6, 15, 44, 30, 47, 38, 37, 54, 60, 13, 51, 39, 48, 50, 3, 10, 21, 52, 24, 57, 28, 27, 32, 17, 20, 40, 12, 29, 59, 45]
        cluster_colors = [cluster_colors[ri] for ri in random_indexes]
        self.cluster_colors = seaborn.color_palette(cluster_colors, as_cmap=True)
    
    ### ------------------------------------ ###
    ### DATA LOADING AND SAVING              ###
    ### ------------------------------------ ###
    
    def import_raw_data(self, input_dir="./", min_counts=500):
        
        # Loading file
        files = [f for f in listdir(input_dir) if f[-3:] == '.h5']

        if not len(files):
            
            print("No HDF5 file found")
            return None

        # Reading the samples data table file or generating one from files prefixes
        try:
            
            samples_table = pd.read_csv('SamplesTable.tsv', index_col=None, header=0, sep="\t")
            
        except:
        
            sample_ids = [f[:f.index('_filtered_feature_bc_matrix.h5')] for f in files]
            samples_table = pd.DataFrame({"SampleID" : sample_ids,
                                          "FileName" : files,
                                          "SampleType" : sample_ids})

        # Reconstructing full matrices and concateneting them
        metadata = {"CellID" : [], "SampleID" : [], "Group" : []}
        for i in range(len(files)):

            new_counts = self.recreate_full_matrix(samples_table.iloc[i,].SampleID, input_dir + '/' + samples_table.iloc[i,].FileName, min_counts)
            metadata["CellID"].extend(new_counts.columns[2:])
            metadata["SampleID"].extend([samples_table.iloc[i,].SampleID for _ in range(len(new_counts.columns) - 2)])
            metadata["Group"].extend([samples_table.iloc[i,].SampleType for _ in range(len(new_counts.columns) - 2)])

            try:
                
                self.raw_counts = pd.concat([self.raw_counts, new_counts.iloc[:, 2:]], axis = 1)
                
            except:
                
                self.raw_counts = new_counts.copy()
        
        self.metadata = pd.DataFrame(metadata)
        
        # Removing features with no counts
        self.raw_counts = self.raw_counts.loc[self.raw_counts.iloc[:, 2:].sum(axis = 1) > 0,]

    ### ------------------------------------ ###
    
    #@staticmethod
    #def load_gtf():
    #    
    #    gtf = pd.DataFrame(data = [line.split('\t') for line in open('genes.gtf').read().split('\n')[5:] if len(line)],
    #                       columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    #    gtf = gtf.loc[gtf.feature == 'exon',]
    #    gtf.attribute = [line.split(';')[0][9 : -1] for line in gtf.attribute]
    #    gtf['feature_length'] = gtf.end.astype(int) - gtf.start.astype(int)
    #    
    #    return gtf
    
    ### ------------------------------------ ###
    
    def load_processed_data(self, file_name, file_dir='./'):
        
        # Reading HDF5 data
        h5_data = h5py.File(file_dir + '/' + file_name, 'r')
        
        # Reconstructing counts matrix
        raw_counts = {}
        raw_counts['GeneID'] = h5_data['features_id'][:].astype(str)
        raw_counts['GeneSymbol'] = h5_data['features_name'][:].astype(str)
        i = 0
        for start, stop in zip(h5_data['cell_range'], h5_data['cell_range'][1:]):
            
            # Extracting cell data
            counts = h5_data['raw_counts'][start : stop]
            indexes = h5_data['cell_features'][start : stop]
            
            # Recreating full counts matrix
            full_counts = np.zeros(len(raw_counts['GeneID']))
            full_counts[indexes] = counts
            
            # Adding data to raw_counts
            cell_name = h5_data['cell_names'][:].astype(str)[i]
            raw_counts[cell_name] = full_counts
            i += 1
        
        self.raw_counts = pd.DataFrame(raw_counts)
        
        # Normalizing and scaling
        self.normalize_counts()
        self.features_mean = h5_data['features_mean'][:]
        self.features_std = h5_data['features_std'][:]
        self.scale_features()
        
        # Reconstructing metadata
        metadata = {}
        metadata['CellID'] = h5_data['cell_names'][:].astype(str).tolist()
        metadata['SampleID'] = h5_data['sample_id'][:].astype(str).tolist()
        metadata['Group'] = h5_data['group'][:].astype(str).tolist()
        
        try:
            
            metadata['Clusters'] = h5_data['clusters'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['Pseudotime'] = h5_data['pseudotime'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['G1PM_Score'] = h5_data['g1pm_score'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['S_Score'] = h5_data['s_score'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['G2M_Score'] = h5_data['g2m_score'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['CellCyclePhase'] = h5_data['cell_cycle_phase'][:].tolist()
            
        except:
            
            pass
        
        try:
            
            metadata['CellCyclePhase'] = h5_data['cell_cycle_phase'][:].astype(str).tolist()
            
        except:
            
            pass
        
        self.metadata = pd.DataFrame(metadata)
        
        # Loading most variable features
        try:
            
            self.most_variable_features = h5_data['most_variable_features'][:].astype(str).tolist()
        
        except:
            
            pass
        
        # Loading PCA and UMAP models
        try:
            
            self.optimal_pca_components = h5_data['optimal_pca_components'][()]
            self.use_variable_features = h5_data['use_variable_features'][()]
            
            # Load PCA model
            self.pca_model = pk.load(open(file_dir + "/pca.pkl", "rb"))
            
            # Load UMAP model(s)
            umap_pk = [file for file in listdir(file_dir) if file[:17] == 'umap_n_neighbors=' and file[-4:] == '.pkl']
            umap_neighbors = [int(file[17:-4]) for file in umap_pk]
            self.umap_models = {n : pk.load(open(file_dir + '/' + model, "rb")) for n,model in zip(umap_neighbors, umap_pk)}
            
            # Reduce dimentions
            self.reduce_dimensions(use_variable_features=self.use_variable_features, pca_components=self.pca_model.n_components_, neighbors=umap_neighbors)
            
        except:
            
            pass
        
        # Loading cluster markers
        try:
            
            cluster_markers = {'Cluster' : h5_data['cluster_markers_cluster'][:].tolist(),
                               'GeneID' : h5_data['cluster_markers_feature_id'][:].astype(str).tolist(),
                               'GeneName' : h5_data['cluster_markers_feature_name'][:].astype(str).tolist(),
                               'ClusterExpression' : h5_data['cluster_markers_cluster_expression'][:].tolist(),
                               'OtherClustersExpression' : h5_data['cluster_markers_other_expression'][:].tolist(),
                               'log2FC' : h5_data['cluster_markers_log2fc'][:].tolist(),
                               'pval' : h5_data['cluster_markers_pval'][:].tolist(),
                               'padj' : h5_data['cluster_markers_padj'][:].tolist()}
            self.cluster_markers = pd.DataFrame(cluster_markers)
        
        except:
            
            pass
        
        # Loading cluster identities
        try:
            
            cluster_identities = {'Cluster' : h5_data['cluster_identities_id'][:].tolist(),
                                  'CellType' : h5_data['cluster_identities_type'][:].astype(str).tolist(),
                                  'CellTypeProb' : h5_data['cluster_identities_type_pval'][:].tolist(),
                                  'CellSubType' : h5_data['cluster_identities_subtype'][:].astype(str).tolist(),
                                  'CellSubTypeProb' : h5_data['cluster_identities_subtype_pval']}
            self.cluster_identities = pd.DataFrame(cluster_identities)
        
        except:
            
            pass
        
        # Reconstructing trajectories
        try:
            
            branches = []
            for br_id in range(h5_data['branches_id'][:].max() + 1):
            
                new_branch = h5_data['branch_vertex'][h5_data['branches_id'][:] == br_id].tolist()
                branches.append(new_branch)
        
            self.branches = branches
            
        except:
            
            pass
        
        # Closing HDF5 file
        h5_data.close()
    
    ### ------------------------------------ ###
    
    def save_class_object(self, object_name="scRNASeq_Analysis"):
        
        if object_name[-3:] != '.h5' or object_name[-5:] != '.hdf5':
            
            object_name += '.h5'
            
        # Raw counts, metadata, optimal pca components, and branches will be stored in a single HDF5 file
        
        # Create sparse matrix for raw features. No need to do this for normalized and scaled features, they'll be created once the object is imported
        cell_range, raw_counts, cell_features = [0], [], []
        for cell_num,cell in enumerate(self.metadata.CellID.to_list()):
            
            # Extracting raw counts for the cell and the index of non-zero features
            cell_counts = self.raw_counts[cell]
            good_features = np.where(cell_counts > 0)[0].tolist()
            cell_counts = cell_counts.iloc[good_features].to_list()
            
            # Storing info
            raw_counts.extend(cell_counts)
            cell_features.extend(good_features)
            cell_range.append(len(raw_counts))
        
        # Storing sample info to HDF5 file
        compressed_data = h5py.File(object_name, 'w')
        
        compressed_data.create_dataset('features_id', data=self.raw_counts['GeneID'].to_list())
        compressed_data.create_dataset('features_name', data=self.raw_counts['GeneSymbol'].to_list())
        
        compressed_data.create_dataset('cell_range', data=cell_range)
        compressed_data.create_dataset('cell_features', data=cell_features)
        compressed_data.create_dataset('raw_counts', data=raw_counts)
        compressed_data.create_dataset('features_mean', data=self.features_mean.tolist())
        compressed_data.create_dataset('features_std', data=self.features_std.tolist())
        
        compressed_data.create_dataset('most_variable_features', data=self.most_variable_features)
        
        compressed_data.create_dataset('optimal_pca_components', data=self.optimal_pca_components)
        compressed_data.create_dataset('use_variable_features', data=self.use_variable_features)
        
        compressed_data.create_dataset('cell_names', data=self.metadata['CellID'].to_list())
        compressed_data.create_dataset('sample_id', data=self.metadata['SampleID'].to_list())
        compressed_data.create_dataset('group', data=self.metadata['Group'].to_list())
        
        if 'Clusters' in self.metadata.columns:
            
            compressed_data.create_dataset('clusters', data=self.metadata['Clusters'].to_list())
            
        if 'Pseudotime' in self.metadata.columns:
            
            compressed_data.create_dataset('pseudotime', data=self.metadata['Pseudotime'].to_list())
        
        if 'G1PM_Score' in self.metadata.columns:
            
            compressed_data.create_dataset('g1pm_score', data=self.metadata['G1PM_Score'].to_list())
        
        if 'S_Score' in self.metadata.columns:
            
            compressed_data.create_dataset('s_score', data=self.metadata['S_Score'].to_list())
        
        if 'G2M_Score' in self.metadata.columns:
            
            compressed_data.create_dataset('g2m_score', data=self.metadata['G2M_Score'].to_list())
        
        if 'CellCyclePhase' in self.metadata.columns:
            
            compressed_data.create_dataset('cell_cycle_phase', data=self.metadata['CellCyclePhase'].to_list())
        
        try:
            
            branches_id = [br_id for br_id,br in enumerate(self.branches) for _ in range(len(br))]
            branch_vertexes = [vertex for br in self.branches for vertex in br]
            compressed_data.create_dataset('branches_id', data=branches_id)
            compressed_data.create_dataset('branch_vertex', data=branch_vertexes)
            
        except:
            
            pass
        
        try:
            
            compressed_data.create_dataset('cluster_identities_id', data=self.cluster_identities.Cluster.to_list())
            compressed_data.create_dataset('cluster_identities_type', data=self.cluster_identities.CellType.to_list())
            compressed_data.create_dataset('cluster_identities_type_pval', data=self.cluster_identities.CellTypeProb.to_list())
            compressed_data.create_dataset('cluster_identities_subtype', data=self.cluster_identities.CellSubType.to_list())
            compressed_data.create_dataset('cluster_identities_subtype_pval', data=self.cluster_identities.CellSubTypeProb.to_list())
            
        except:
            
            pass

        try:
            
            compressed_data.create_dataset('cluster_markers_cluster', data=self.cluster_markers.Cluster.to_list())
            compressed_data.create_dataset('cluster_markers_feature_id', data=self.cluster_markers.GeneID.to_list())
            compressed_data.create_dataset('cluster_markers_feature_name', data=self.cluster_markers.GeneName.to_list())
            compressed_data.create_dataset('cluster_markers_cluster_expression', data=self.cluster_markers.ClusterExpression.to_list())
            compressed_data.create_dataset('cluster_markers_other_expression', data=self.cluster_markers.OtherClustersExpression.to_list())
            compressed_data.create_dataset('cluster_markers_log2fc', data=self.cluster_markers.log2FC.to_list())
            compressed_data.create_dataset('cluster_markers_pval', data=self.cluster_markers.pval.to_list())
            compressed_data.create_dataset('cluster_markers_padj', data=self.cluster_markers.padj.to_list())
            
        except:
            
            pass
        
        compressed_data.close()
        
        # Saving PCA and UMAP models as Pickle files
        try:
            
            pk.dump(self.pca_model, open("pca.pkl", "wb"))
            
        except:
            
            pass
        
        try:
            
            for n,umap_n in self.umap_models.items():
                
                pk.dump(umap_n, open("umap_n_neighbors={}.pkl".format(n), "wb"))
                
        except:
            
            pass

    ### ------------------------------------ ###

    @staticmethod
    def recreate_full_matrix(file_id, file, min_counts):
    
        # Reading HDF5 data
        h5_data = h5py.File(file, 'r')

        # Initializing empty matrix
        feature_types = np.array(h5_data['matrix']['features']['feature_type']).astype(str)
        gene_ids = np.array(h5_data['matrix']['features']['id']).astype(str)
        gene_names = np.array(h5_data['matrix']['features']['name']).astype(str)
        counts_matrix = {'GeneID' : gene_ids, 'GeneSymbol' : gene_names}

        # Adding cells (if their feature counts >= min_counts)
        cell_count = 1

        for start, stop in zip(h5_data['matrix']['indptr'], h5_data['matrix']['indptr'][1:]):

            # Extracting cell data
            counts = h5_data['matrix']['data'][start : stop]
            indexes = h5_data['matrix']['indices'][start : stop]

            if counts.sum() < min_counts:
            
                continue

            # Recreating full counts matrix
            full_counts = np.zeros(len(gene_ids))
            full_counts[indexes] = counts

            # Adding data to counts_matrix
            leading_zeroes = ''.join(['0' for limit in [10, 100, 1000, 10000] if cell_count < limit])
            cell_name = '{}_{}{}'.format(file_id, leading_zeroes, cell_count)
            counts_matrix[cell_name] = full_counts

            cell_count += 1

        # Closing HDF5 file
        h5_data.close()

        # Converting to dataframe
        counts_matrix = pd.DataFrame(counts_matrix)

        # Demultiplexing
        if 'Antibody Capture' in feature_types:
            
            # Extracting hashtags counts
            hashtag_num = np.unique(gene_ids[feature_types == 'Antibody Capture'])
            hashtag_ids = np.unique(gene_names[feature_types == 'Antibody Capture'])
            
            hashtag_counts = counts_matrix.loc[counts_matrix.GeneID.isin(hashtag_num),]
            
            counts_matrix = counts_matrix.loc[~counts_matrix.GeneID.isin(hashtag_num)]
            
            # Normalizing counts (centered log ratio)
            log_geometric_means = [(1 / len(h_counts)) * np.log(h_counts[h_counts != 0]).sum() for h_counts in hashtag_counts.iloc[:, 2:].to_numpy()]
            
            for i in range(hashtag_counts.shape[0]):
                
                normalized_hash_counts = hashtag_counts.iloc[i, 2:].to_numpy()
                normalized_hash_counts[normalized_hash_counts != 0] = np.log(normalized_hash_counts[normalized_hash_counts != 0.].tolist()) - log_geometric_means[i]
                
                hashtag_counts.iloc[i, 2:] = normalized_hash_counts.copy()
            
            # Sorting geometric means, then if from mean i to i+1 there's a drop >= 75%, discard i+i and all that follows
            sorted_geometric_means = np.sort(np.exp(log_geometric_means))[::-1]
            threshold = [sgm1 - 1 for sgm1,sgm2 in zip(sorted_geometric_means[:-1], sorted_geometric_means[1:]) if sgm2 / sgm1 <= 0.25]
            threshold = threshold[-1] if len(threshold) else min(sorted_geometric_means) - 1
            
            good_hashtag_ids = np.array([h_id for h_id,gm in zip(hashtag_ids, np.exp(log_geometric_means)) if gm > threshold])
            
            hashtag_counts = hashtag_counts.loc[hashtag_counts.GeneSymbol.isin(good_hashtag_ids)]
            
            # Initial clustering
            n_clusters = len(good_hashtag_ids) + 1
            cluster_labels = KMeans(n_clusters).fit_predict(hashtag_counts.iloc[:, 2:].T)
            
            # Computing hashtag levels threshold
            classification = {}
            for hashtag in good_hashtag_ids:
                
                # Computing clusters average levels
                clusters_average = []
                
                for cl in range(n_clusters):

                    cl_avg = hashtag_counts.loc[hashtag_counts.GeneSymbol == hashtag, [False, False] + (cluster_labels == cl).tolist()].to_numpy().mean()
                    clusters_average.append(cl_avg)

                # Using the cluster with lowest average levels for fitting a negative binomial distribution, then calculating the quantile at q = probability of failure = 0.99
                cl = np.argmin(clusters_average)
                
                cell_subset = hashtag_counts.loc[hashtag_counts.GeneSymbol == hashtag, [False, False] + (cluster_labels == cl).tolist()].to_numpy()[0]
                #threshold = np.percentile(cell_subset, 99) # A simpler way is to set the threshold as the 99th percentile
                
                cell_subset[cell_subset < 0] = 0 # Only positive values are considered by a neagtive binomial distribution
                
                log_mu, alpha = NegativeBinomial(cell_subset, np.ones_like(cell_subset)).fit(skip_hessian=True, disp=0).params
                mu = np.exp(log_mu)
                p = 1 / (1 + mu * alpha)
                n = mu * p / (1 - p)
                threshold = nbinom.ppf(0.99, n, p, loc=mu)
                
                classification[hashtag] = (hashtag_counts.loc[hashtag_counts.GeneSymbol == hashtag,].iloc[0, 2:] > threshold).to_list()
            
            classification = pd.DataFrame(classification)

            # Reporting results of demultiplexing
            negatives, singlets, multiplets = sum(classification.sum(axis=1) == 0), sum(classification.sum(axis=1) == 1), sum(classification.sum(axis=1) > 1)
            tot = negatives + singlets + multiplets
            print(f'Demultiplexing result:\n{singlets} ({round(singlets / tot, 3)}%) cells were singlets.\n{negatives + multiplets} ({round((negatives + multiplets) / tot, 3)}%) cells were either negative or multiplets')

            # Removing non-singlets
            keep_cells = (classification.sum(axis=1) == 1).to_numpy()
            classification = classification.loc[keep_cells,]
            counts_matrix = counts_matrix.loc[:, [True, True] + keep_cells.tolist()]
            
            # Renaming cells to add hashtag classification
            new_names = ['GeneID', 'GeneSymbol']
            for cell_count in range(classification.shape[0]):
                
                leading_zeroes = ''.join(['0' for limit in [10, 100, 1000, 10000] if cell_count < limit])
                identity = good_hashtag_ids[np.argmax(classification.iloc[cell_count,])] if sum(classification.iloc[cell_count,]) == 1 else ''
            
                cell_name = '{}_{}_{}{}'.format(file_id, identity, leading_zeroes, cell_count)
                new_names.append(cell_name)
            
            counts_matrix.columns = new_names

        return counts_matrix
    
    ### ------------------------------------ ###
    ### DATA PREPROCESSING                   ###
    ### ------------------------------------ ###
    
    def find_most_variable_features(self, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000):
        
        # Graph based, similar to Seurat's FindVariableGenes
        
        # Calculating mean and log(variance/mean) (X and Y respectively) for each feature
        var_plot = pd.DataFrame({'GeneID' : self.normalized_counts.GeneID,
                                 'X' : np.exp(self.normalized_counts.iloc[:, 2:].mean(axis = 1)),
                                 'Y' : np.log(self.normalized_counts.iloc[:, 2:].var(axis = 1) / (self.normalized_counts.iloc[:, 2:].mean(axis = 1) + 0.00001))})
        min_x, max_x = max(x_low_cutoff, var_plot.X.min()), min(x_high_cutof, var_plot.X.max())
        bin_width = (max_x - min_x) / nbins

        # Binning and computing log(variance/mean) z-scores
        most_variable_features = {'GeneID' : [], 'ZScores' : []}

        for x_start in np.arange(min_x, max_x, bin_width):
            
            x_stop = x_start + bin_width
            
            binned = var_plot.loc[(var_plot.X > x_start) & (var_plot.X <= x_stop),]
            
            var_z_scores = (binned.Y - binned.Y.mean()) / binned.Y.std()

            most_variable_features['GeneID'].extend(list(binned.GeneID))
            most_variable_features['ZScores'].extend(list(var_z_scores))

        most_variable_features = pd.DataFrame(most_variable_features)
        most_variable_features.sort_values(by = 'ZScores', axis = 0, ascending = False, inplace = True)

        # Removing features with too low/high z_score, then selecting top ones
        self.most_variable_features = list(most_variable_features.loc[(most_variable_features.ZScores >= y_low_cutoff) &
                                                                      (most_variable_features.ZScores <= y_high_cutoff),].iloc[:max_features,].GeneID)
    
    ### ------------------------------------ ###
    
    def normalize_counts(self, scale_factor=10000):
        
        # Copying raw_counts matrix
        self.normalized_counts = self.raw_counts.copy()
        
        # LogNormalization, similar to Seurat's NormalizeData
        scale_factor = 10000
        counts_per_cell = np.array(self.normalized_counts.iloc[:, 2:].sum(axis = 0))
        self.normalized_counts.iloc[:, 2:] = np.log1p(self.normalized_counts.iloc[:, 2:].div(counts_per_cell / scale_factor, axis = 1))
    
    ### ------------------------------------ ###

    def scale_features(self):
        
        # Copying normalized_counts matrix
        self.scaled_features = self.normalized_counts.copy()
        
        # Scaling features
        try: # Precalculated means and stds
            
            self.scaled_features.iloc[:, 2:] = self.scaled_features.iloc[:, 2:].subtract(self.features_mean, axis="index").div(self.features_std, axis="index")
        
        except:
        
            self.features_mean = self.normalized_counts.iloc[:, 2:].mean(axis = 1).to_numpy()
            self.features_std = self.normalized_counts.iloc[:, 2:].std(axis = 1).to_numpy()
            self.scaled_features.iloc[:, 2:] = self.scaled_features.iloc[:, 2:].subtract(self.features_mean, axis="index").div(self.features_std, axis="index")

    ### ------------------------------------ ###
    ### DIMENSIONALITY REDUCTION             ###
    ### ------------------------------------ ###

    def reduce_dimensions(self, use_variable_features=True, pca_components=50, neighbors=[5, 30, 50]):
        
        # Using either only the most_variable_features or all of them
        self.use_variable_features = use_variable_features
        if self.use_variable_features:
            
            pca_input = self.scaled_features.loc[self.scaled_features.GeneID.isin(self.most_variable_features),].iloc[:, 2:].to_numpy().T
            
        else:
            
            pca_input = self.scaled_features.iloc[:, 2:].to_numpy().T
        
        try:
            
            self.pca = self.pca_model.transform(pca_input)
            
        except:
            
            self.pca_model = PCA(pca_components)
            self.pca = self.pca_model.fit_transform(pca_input)
        
            # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
            x0, x1 = 0, min(len(self.pca_model.explained_variance_), pca_components - 1)
            y0, y1 = self.pca_model.explained_variance_[x0], self.pca_model.explained_variance_[x1]
            gradient = (y1 - y0) / (x1 - x0)
            intercept = y0
            difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(self.pca_model.explained_variance_[:pca_components])]
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
        
        # Running UMAP with X neighbors
        try:
            
            self.umaps = {}
            for n in neighbors:
                
                counts_umap = self.umap_models[n].transform(self.pca[:, :self.optimal_pca_components])
                counts_umap = pd.DataFrame(data = counts_umap, index = self.raw_counts.columns[2:], columns = ['UMAP_1', 'UMAP_2'])
                
                self.umaps[n] = counts_umap
        
        except:
            
            self.umap_models = {}
            self.umaps = {}
            for n in neighbors:
            
                umap_n = umap.UMAP(n_components=2, n_neighbors=n, random_state=42)
                counts_umap = umap_n.fit_transform(self.pca[:, :self.optimal_pca_components])
                counts_umap = pd.DataFrame(data = counts_umap, index = self.raw_counts.columns[2:], columns = ['UMAP_1', 'UMAP_2'])
                
                self.umap_models[n] = umap_n
                self.umaps[n] = counts_umap
    
    ### ------------------------------------ ###
    ### CLUSTERING                           ###
    ### ------------------------------------ ###

    def cluster_cells(self, n_neighbors=10):
        
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
        
        # Updating metadata
        self.metadata["Clusters"] = clusters
     
    ### ------------------------------------ ###
    
    def compare_populations(self, pop1, pop2, min_fc=0.25, min_pct=0.1):
        
        # Extracting average exression values and log2Fc
        pop1_expression = self.normalized_counts.loc[:, pop1].mean(axis = 1)
        pop2_expression = self.normalized_counts.loc[:, pop2].mean(axis = 1)
        log2_fc = np.log2(pop1_expression / (pop2_expression + 0.001) + 0.001)
        
        # Calculating the number of cells in which a certain feature is detected (i.e. normalized_count > 0)
        pct_pop1 = (self.normalized_counts.loc[:, pop1] > 0).sum(axis=1) / len(pop1)
        pct_pop2 = (self.normalized_counts.loc[:, pop2] > 0).sum(axis=1) / len(pop2)
        
        # Filtering features to test by min_fc and min_pct
        features_to_test = (abs(log2_fc) > min_fc) & ((pct_pop1 > min_pct) | (pct_pop2 > min_pct))
        
        # Significance testing and BH correction
        pvals = mannwhitneyu(self.normalized_counts.loc[features_to_test, pop1], self.normalized_counts.loc[features_to_test, pop2], alternative="two-sided", axis=1).pvalue
        _, padjs = fdrcorrection(pvals, alpha=0.05, is_sorted=False)
        
        # Filtering
        significant_genes_id = list(self.normalized_counts.GeneID[features_to_test][padjs < 0.05])
        significant_genes_name = list(self.normalized_counts.GeneSymbol[features_to_test][padjs < 0.05])
        log2_fc = log2_fc[features_to_test][padjs < 0.05]
        pop1_expression = pop1_expression[features_to_test][padjs < 0.05]
        pop2_expression = pop2_expression[features_to_test][padjs < 0.05]
        pvals = pvals[padjs < 0.05]
        padjs = padjs[padjs < 0.05]
        
        return significant_genes_id, significant_genes_name, log2_fc, pop1_expression, pop2_expression, pvals, padjs

    ### ------------------------------------ ###
    
    def find_clusters_markers(self, min_fc=0.25, min_pct=0.1):
        
        cluster_markers = {"Cluster" : [], "GeneID" : [], "GeneName" : [], "ClusterExpression" : [], "OtherClustersExpression" : [], "log2FC" : [], "pval" : [], "padj" : []}
        
        clusters = self.metadata.Clusters.unique()
        clusters.sort()
        
        for cl in clusters:
            
            print("Finding markers of cluster {}".format(cl))
            
            cluster_cells = list(self.metadata.loc[self.metadata.Clusters == cl, "CellID"])
            other_cells = list(self.metadata.loc[self.metadata.Clusters != cl, "CellID"])
            
            significant_genes_id, significant_genes_name, log2_fc, pop1_expression, pop2_expression, pvals, padjs = self.compare_populations(cluster_cells, other_cells, min_fc, min_pct)
            
            cluster_markers["Cluster"].extend([cl for _ in range(len(significant_genes_id))])
            cluster_markers["GeneID"].extend(significant_genes_id)
            cluster_markers["GeneName"].extend(significant_genes_name)
            cluster_markers["ClusterExpression"].extend(pop1_expression)
            cluster_markers["OtherClustersExpression"].extend(pop2_expression)
            cluster_markers["log2FC"].extend(log2_fc)
            cluster_markers["pval"].extend(pvals)
            cluster_markers["padj"].extend(padjs)
        
        self.cluster_markers = pd.DataFrame(cluster_markers)
    
    ### ------------------------------------ ###
    
    def score_cluster_identity(self, ref_expression, ref_samples, clusters=[]):
        
        # Similar to CIPR (see https://github.com/atakanekiz/CIPR-Package), but in Python and providing softmax normalized scores

        # Listing clusters used
        if not len(clusters):
        
            clusters = self.metadata.Clusters.unique()
            clusters.sort()

        # Listing possible cell types and subtypes
        cell_classes = ref_samples.reference_cell_type
        cell_subclasses = ref_samples.long_name
        
        # Init cluster identity matrix
        cluster_identities = {'Cluster' : [],
                              'CellType' : [],
                              'CellTypeProb' : [],
                              'CellSubType' : [],
                              'CellSubTypeProb' : []}
        
        # Assessing identity probabilities for each cluster
        for cl in clusters:
            
            # Subset markers
            cl_markers = self.cluster_markers.loc[(self.cluster_markers.Cluster == cl) &
                                                  (self.cluster_markers.padj < 0.05) &
                                                  (self.cluster_markers.GeneName.str.lower().isin(ref_expression.GeneName.str.lower())),]
            cl_markers = cl_markers.loc[cl_markers.GeneName.str.lower().duplicated() == False,]
            cl_markers = cl_markers.iloc[np.argsort(cl_markers.GeneName.str.lower()),]
            
            # Testing cell_types
            probs = []
            for col in ref_expression.columns[1:]:
                
                ref_logfc = ref_expression.loc[ref_expression.GeneName.str.lower().isin(cl_markers.GeneName.str.lower()), ["GeneName", col]]
                ref_logfc = ref_logfc.iloc[np.argsort(ref_logfc.GeneName.str.lower()),]
                
                score = (np.array(cl_markers.log2FC) * np.array(ref_logfc.loc[:,col])).sum()
                
                probs.append(score)
            
            # Scaling scores to -100:100 range
            probs = np.array(probs)
            probs = ((probs - probs.min()) / (probs.max() - probs.min())) * 200 - 100
        
            # Softmax
            probs = np.exp(probs) / sum(np.exp(probs))
            
            # Finding most probable class
            best_class = ""
            best_class_prob = 0
            for c_class in cell_classes:
                
                class_prob = probs[ref_samples.reference_cell_type == c_class].sum()
                
                if class_prob > best_class_prob:
                    
                    best_class = c_class
                    best_class_prob = class_prob

            # Finding most probable subclass
            subclass_probs = probs[ref_samples.reference_cell_type == best_class]
            best_subclass_prob = subclass_probs.max()
            # FIX THIS
            best_subclass = list(cell_subclasses[ref_samples.reference_cell_type == best_class][subclass_probs == best_subclass_prob])[0]
            
            # Adding cluster info
            cluster_identities['Cluster'].append(cl)
            cluster_identities['CellType'].append(best_class)
            cluster_identities['CellTypeProb'].append(best_class_prob)
            cluster_identities['CellSubType'].append(best_subclass)
            cluster_identities['CellSubTypeProb'].append(best_subclass_prob)
        
        self.cluster_identities = pd.DataFrame(cluster_identities)
    
    ### ------------------------------------ ###
    ### TRAJECTORY ANALYSIS                  ###
    ### ------------------------------------ ###
    
    def create_mst(self, clusters=[]):
        
        # Kruskal Minimal Spanning Tree algorithm
        
        # Finding clusters centroids (graph vertices)
        if not len(clusters):
        
            clusters = list(self.metadata.Clusters.unique())
            clusters.sort()
        
        centroids = np.array([self.pca[self.metadata.Clusters == cl, :self.optimal_pca_components].mean(axis = 0) for cl in clusters])
        
        # Calculating pairwise distances between cluster centroids (graph edges)
        distances = []
        for c1_index,c1 in enumerate(centroids):
            
            for c2_index,c2 in enumerate(centroids[c1_index+1:]):
                
                new_distance = (((c1 - c2)**2).sum())**0.5
                distances.append([c1_index, c1_index + c2_index + 1, new_distance])
        
        distances = np.array(distances)
        distances = distances[distances[:, 2].argsort(),]
        
        # Init parents and ranks
        n_vertices = len(clusters)
        parents = list(range(n_vertices))
        ranks = [0 for _ in range(n_vertices)]
        
        # Generating MST
        good_edges = []
        edge_count =  0
        current_edge = -1
        
        while edge_count < n_vertices - 1:
            
            current_edge += 1
    
            # Extracting vertices
            vert_1, vert_2 = distances[current_edge, :2]
            
            # Finding parent of vertex 1
            parent_1 = int(vert_1)
            while parents[parent_1] != parent_1:
      
                parent_1 = parents[parent_1]
            
            # Finding parent of vertex 2
            parent_2 = int(vert_2)
            while parents[parent_2] != parent_2:
      
                parent_2 = parents[parent_2]
            
            if parent_1 != parent_2:
      
                good_edges.append(current_edge)
                edge_count += 1
      
                # Union
                if ranks[parent_1] > ranks[parent_2]:
        
                    parents[parent_2] = parent_1
        
                elif ranks[parent_2] > ranks[parent_1]:
        
                    parents[parent_1] = parent_2
        
                else:
        
                    parents[parent_2] = parent_1
                    ranks[parent_1] += 1
        
        # Updating distance matrix
        distances = distances[good_edges,]
        
        # Finding the parent of all
        updated_parents = []
        for i,p in enumerate(parents):
            
            parent = p
            
            while parents[parent] != parent:
                
                parent = parents[parent]
                
            updated_parents.append(parent)
        
        odin = list(set(updated_parents))[0] # Odin the allfather!
        
        # Reconstructing branches by moving out from odin vertex
        # Branches are initialized using the edges directly connected to odin
        used_edges = np.where(distances[:,:2] == odin)[0].tolist()
        branches = [[odin, int(d[0])] if int(d[1]) == odin else [odin, int(d[1])] for d in distances[used_edges] if odin in d]
        
        while len(used_edges) != len(distances):
            
            old_branches = [br[:] for br in branches]
            
            for br_index,br in enumerate(old_branches):
                
                # Finding edges where br[-1] is a vertex
                rows, cols = np.where(distances[:,:2] == br[-1])
                cols = [1 if c == 0 else 0 for c in cols]
                toggle = True
                
                for r,c in zip(rows, cols):
                    
                    if r not in used_edges:
                        
                        if toggle: # If more than one edge is found, extend the branch once, then create new branches
                        
                            branches[br_index].append(int(distances[r,c]))
                            toggle = False
                            
                        else:
                            
                            branches.append(br + [int(distances[r,c])])
                            
                        used_edges.append(r)
        
        # Converting branches from cluster index to cluster name
        for br_n,br in enumerate(branches):
            
            branches[br_n] = [clusters[index] for index in br]
                    
        self.branches = branches

    ### ------------------------------------ ###
    
    def generate_pseudotime(self):
        
        # Listing clusters used
        clusters = list({vertex for br in self.branches for vertex in br})
        clusters.sort()
        
        # Calculating cluster centroids
        centroids = np.array([self.pca[self.metadata.Clusters == cl, :self.optimal_pca_components].mean(axis = 0) for cl in clusters])
        
        # Calculating pseudotime for each cell
        pseudotimes = []
        for c in range(len(self.metadata)):
            
            coords = self.pca[c, :self.optimal_pca_components]
            cluster = list(self.metadata.Clusters)[c]
            
            branch = [br for br in self.branches if cluster in br]
            
            if not len(branch) or cluster not in clusters: # Cluster is not connected
                
                pseudotimes.append(0)
                continue
            
            branch = branch[0]
            
            vertex_index = branch.index(cluster)
            
            if vertex_index == 0: # Cell belongs to root vertex. Pseudotime = projection on first edge of branch
                
                # Projecting cell to first edge of branch
                branch_centroids = np.array([centroids[clusters.index(cl)] for cl in branch[:2]])
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
                
                # Computing pseudotime
                cell_pseudotime = ((projection**2).sum())**0.5
                
            elif vertex_index == len(branch) - 1: # Cell belongs to end vertex. Pseudotime = projection on last edge of branch
                
                # Projecting cell to last edge of branch
                branch_centroids = np.array([centroids[clusters.index(cl)] for cl in branch[-2:]])
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
          
                # Calculating previous branch length
                previous_centroids = np.array([centroids[clusters.index(cl)] for cl in branch[:vertex_index]])
                added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
          
                # Computing pseudotime
                cell_pseudotime = ((projection**2).sum())**0.5 + added_distance
                
            else: # Cell belongs to a end/middle vertex. Pseudotime = projection on edge of branch + previous branches' length
                
                # Finding closest cluster centroid to the cell that the cell doesn't belong to
                near_centroids = np.array([centroids[clusters.index(cl)] for cl in branch[vertex_index - 1 : vertex_index + 3]])
                distance_1 = (((near_centroids[0] - coords)**2).sum())**0.5
                distance_2 = (((near_centroids[-1] - coords)**2).sum())**0.5
                
                # Projecting cell to nearest edge of branch
                branch_centroids = near_centroids[:2] if distance_1 < distance_2 else near_centroids[-2:]
                vertex_index = vertex_index if distance_1 < distance_2 else vertex_index + 1
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
          
                # Calculating previous branch length
                previous_centroids = np.array([centroids[clusters.index(cl)] for cl in branch[:vertex_index]])
                added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
          
                # Computing pseudotime
                cell_pseudotime = ((projection**2).sum())**0.5 + added_distance
                
            pseudotimes.append(cell_pseudotime)
        
        # Normalize pseudotime to 0-100 and update metadata
        self.metadata["Pseudotime"] = (np.array(pseudotimes) - min(pseudotimes)) * 100 / (max(pseudotimes) - min(pseudotimes))
    
    ### ------------------------------------ ###
    ### CELL CYCLE ANALYSIS                  ###
    ### ------------------------------------ ###

    def score_cell_cycle(self, g1pm_features, s_features, g2m_features, gene_pool=[], bins=50, ctrl_genes_num=50):
        
        phases = ["G1PM", "S", "G2M"]
        
        # Scoring gene sets
        toggle = 1
        cell_cycle_scores = []
        for set_name,gene_set in zip(phases, [g1pm_features, s_features, g2m_features]):
            
            set_score = self.score_gene_set(gene_set, gene_pool, bins, ctrl_genes_num)

            if set_score is None:
                
                toggle = 0
                bad_gene_set = set_name
                break
            
            else:
                
                cell_cycle_scores.append(set_score)
        
        if toggle == 0:
            
            print("Error: {} list is bad".format(bad_gene_set))
            return None
        
        else:
            
            cell_cycle_scores = pd.concat(cell_cycle_scores, axis=1)
            cell_cycle_scores.columns = ["{}_Score".format(p) for p in phases]

            cell_cycle_scores["CellCyclePhase"] = [phases[x] if cell_cycle_scores.iloc[y, x] > 0 else "G1" for y,x in enumerate(np.argmax(cell_cycle_scores.to_numpy(), axis=1))]
            
            self.metadata = self.metadata.merge(cell_cycle_scores, left_on='CellID', right_on=cell_cycle_scores.index)

    ### ------------------------------------ ###

    def score_gene_set(self, gene_set, gene_pool=[], bins=50, ctrl_genes_num=50):

        # Subsetting gene_set and gene_list for detected genes
        gene_set = [gene for gene in gene_set if gene in self.scaled_features.GeneID.tolist() or gene in self.scaled_features.GeneSymbol.tolist()]
        gene_pool = [gene for gene in gene_pool if gene in self.scaled_features.GeneID.tolist() or gene in self.scaled_features.GeneSymbol.tolist()]

        # Making sure that gene_set and gene_pool have len > 0
        if not len(gene_set):
            
            print("Error: no gene_list provided")
            return None
        
        if not len(gene_pool):
            
            gene_pool = list(self.scaled_features.GeneID)
        
        # Mean of each gene in the gene_pool across all cells
        gene_means = self.scaled_features.loc[(self.scaled_features.GeneID.isin(gene_pool)) | (self.scaled_features.GeneSymbol.isin(gene_pool))].iloc[:, 2:].mean(axis=1)
        gene_means.index = gene_pool
        
        # Rank genes based on binned expression level, then for each bin of the genes in the gene_set, pick ctrl_genes random genes for genes with matched binned expression
        bin_size = len(gene_means) // bins
        expression_order = np.argsort(gene_means)
        
        set_indexes = []
        for gene in gene_set:
            
            try:
                index = expression_order.iloc[self.scaled_features.GeneID.tolist().index(gene)]
            except:
                index = expression_order.iloc[self.scaled_features.GeneSymbol.tolist().index(gene)]
            
            set_indexes.append(index)
        
        ctrl_set = []
        for index in set_indexes:
            
            random_pool = expression_order.index[(expression_order >= (index - bin_size // 2)) &
                                                 (expression_order <= (index + bin_size // 2)) &
                                                 ~(expression_order.isin(set_indexes))].tolist()
            np.random.shuffle(random_pool)
            ctrl_set.extend(random_pool[:ctrl_genes_num])
        ctrl_set = list(set(ctrl_set)) # Removing duplicates
        
        # Computing the mean of gene_set and ctrl_set genes for each cell
        set_means = self.scaled_features.loc[(self.scaled_features.GeneID.isin(gene_set)) | (self.scaled_features.GeneSymbol.isin(gene_set))].iloc[:, 2:].mean(axis=0)
        ctrl_means = self.scaled_features.loc[(self.scaled_features.GeneID.isin(ctrl_set)) | (self.scaled_features.GeneSymbol.isin(ctrl_set))].iloc[:, 2:].mean(axis=0)
        
        set_score = set_means - ctrl_means
        
        return set_score
    
    ### ------------------------------------ ###
    ### PLOTTING                             ###
    ### ------------------------------------ ###

    def plot_cell_cycle(self, cells=[], samples=[], groups=[], clusters=[], dot_size=1.5):
        
        # Creating a filter for cells of interest
        if not len(cells):
        
            cells = self.metadata.CellID
        
        if not len(samples):
        
            samples = self.metadata.CellID.unique()
        
        if not len(groups):
        
            groups = self.metadata.Group.unique()

        if not len(clusters):
        
            clusters = np.array([True for _ in range(len(self.metadata))])
        
        cell_filter = list(self.metadata.CellID.isin(cells) & self.metadata.CellID.isin(samples) & self.metadata.Group.isin(groups) & clusters)
        
        if "CellCyclePhase" not in self.metadata.columns:
            
            print("Please run score_cell_cycle() before plotting")
            
        else:
            
            for i, (n, coords) in enumerate(self.umaps.items()):

                # Adding expression data, then sorting by smallest value
                plot_data = coords.copy()
                plot_data['CellCyclePhase'] = list(self.metadata["CellCyclePhase"])
                
                # Subsetting cells
                plot_data = plot_data.loc[cell_filter,]
                
                # Plotting
                plt.figure(figsize=(5, 5))
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='CellCyclePhase', hue_order=['G1PM', 'G1', 'S', 'G2M'], palette='Set1', marker='.', s=dot_size, linewidth=0)
                legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Cell Cycle Phase')
                plt.xlabel('UMAP 1')
                plt.ylabel('UMAP 2')
                plt.savefig("Cell_Cycle_Phase_UMAP-{}.png".format(n), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
                plt.close()

    ### ------------------------------------ ###

    def plot_clusters(self, cells=[], samples=[], groups=[], clusters=[], dot_size=1.5):
        
        # Creating a filter for cells of interest
        if not len(cells):
        
            cells = self.metadata.CellID
        
        if not len(samples):
        
            samples = self.metadata.CellID.unique()
        
        if not len(groups):
        
            groups = self.metadata.Group.unique()

        if not len(clusters):
        
            clusters = list(self.metadata.Clusters.unique())
            clusters.sort()
        
        cell_filter = list(self.metadata.CellID.isin(cells) & self.metadata.CellID.isin(samples) & self.metadata.Group.isin(groups) & self.metadata.Clusters.isin(clusters))

        # Subsetting cluster color palette
        color_palette = [self.cluster_colors[cl] for cl in clusters]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(clusters) / 16) + 1

        for i, (n, coords) in enumerate(self.umaps.items()):

            # Adding expression data, then sorting by smallest value
            plot_data = coords.copy()
            plot_data['Clusters'] = list(self.metadata.Clusters)
        
            # Subsetting cells
            plot_data = plot_data.loc[cell_filter,]
        
            # Finding cluster centers for writing text
            centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
        
            # Plotting
            plt.figure(figsize=(5, 5))
            seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Clusters', palette=color_palette, marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Clusters', ncol=legend_cols)
            for c_num,cl in enumerate(centers):
                x, y = cl
                plt.text(x, y, str(clusters[c_num]), horizontalalignment='center', size='small', color='black', weight='semibold')
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig("Clusters_UMAP-{}.png".format(n), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
            plt.close()
    
    ### ------------------------------------ ###
    
    def plot_gene_expression(self, target_gene, cells=[], samples=[], groups=[], clusters=[], dot_size=1.5):
        
        # Creating a filter for cells of interest
        if not len(cells):
        
            cells = self.metadata.CellID
        
        if not len(samples):
        
            samples = self.metadata.CellID.unique()
        
        if not len(groups):
        
            groups = self.metadata.Group.unique()

        if not len(clusters):
        
            clusters = np.array([True for _ in range(len(self.metadata))])
        
        cell_filter = list(self.metadata.CellID.isin(cells) & self.metadata.CellID.isin(samples) & self.metadata.Group.isin(groups) & clusters)
        
        # Extracting expression data for target_gene
        if target_gene not in list(self.raw_counts.GeneID) and target_gene not in list(self.raw_counts.GeneSymbol):
            
            toggle = False
        
        else:
            
            toggle = True
            
            try:
                
                gene_index = list(self.raw_counts.GeneID).index(target_gene)
                
            except:
                
                gene_index = list(self.raw_counts.GeneSymbol).index(target_gene)
                
            expression = self.normalized_counts.iloc[gene_index, 2:].to_numpy()
        
        if toggle:
            
            for i, (n, coords) in enumerate(self.umaps.items()):

                # Adding expression data
                plot_data = coords.copy()
                plot_data['Expression'] = expression
                
                # Subsetting cells
                plot_data = plot_data.loc[cell_filter,]
                
                # Plotting
                plt.figure(figsize=(5, 5))
                if plot_data.Expression.min() != plot_data.Expression.max():
                    seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Expression', palette='viridis', marker='.', s=dot_size, linewidth=0)
                else:
                    seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Expression', hue_norm=(0, 100), palette='viridis', marker='.', s=dot_size, linewidth=0)
                legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='{} Expression'.format(target_gene))
                plt.xlabel('UMAP 1')
                plt.ylabel('UMAP 2')
                plt.savefig("{}_Expression_UMAP-{}.png".format(target_gene, n), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
                plt.close()
        
        else:
            
            print("Target gene {}, not found".format(target_gene))
    
    ### ------------------------------------ ###

    def plot_gene_set_score(self, set_score, gene_set_name='', cells=[], samples=[], groups=[], clusters=[], dot_size=1.5):
        
        # Check gene set name
        if gene_set_name == '':

            gene_set_name = 'Gene_Set'

        # Creating a filter for cells of interest
        if not len(cells):
        
            cells = self.metadata.CellID
        
        if not len(samples):
        
            samples = self.metadata.CellID.unique()
        
        if not len(groups):
        
            groups = self.metadata.Group.unique()

        if not len(clusters):
        
            clusters = np.array([True for _ in range(len(self.metadata))])
        
        cell_filter = list(self.metadata.CellID.isin(cells) & self.metadata.CellID.isin(samples) & self.metadata.Group.isin(groups) & clusters)
        
        for i, (n, coords) in enumerate(self.umaps.items()):

            # Adding expression data, then sorting by smallest value
            plot_data = coords.copy()
            plot_data['Score'] = list(set_score)
            
            # Subsetting cells
            plot_data = plot_data.loc[cell_filter,]
            
            # Plotting
            plt.figure(figsize=(5, 5))
            if plot_data.Score.min() != plot_data.Score.max():
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Score', palette='viridis', marker='.', s=dot_size, linewidth=0)
            else:
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Score', hue_norm=(0, 100), palette='viridis', marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Gene Set Expression')
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig("{}_Expression_UMAP-{}.png".format(gene_set_name, n), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
            plt.close()

    ### ------------------------------------ ###

    def plot_trajectories(self, cells=[], samples=[], groups=[], pseudotime=False, dot_size=1.5):
        
        # Creating a filter for cells of interest
        if not len(cells):
        
            cells = self.metadata.CellID
        
        if not len(samples):
        
            samples = self.metadata.CellID.unique()
        
        if not len(groups):
        
            groups = self.metadata.Group.unique()

        clusters = list({vertex for br in self.branches for vertex in br})
        clusters.sort()
        
        cell_filter = list(self.metadata.CellID.isin(cells) & self.metadata.CellID.isin(samples) & self.metadata.Group.isin(groups) & self.metadata.Clusters.isin(clusters))
        
        # Subsetting cluster color palette
        color_palette = [self.cluster_colors[cl] for cl in clusters]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(clusters) / 16) + 1

        for i, (n, coords) in enumerate(self.umaps.items()):

            # Adding expression data, then sorting by smallest value
            plot_data = coords.copy()
            plot_data['Clusters'] = list(self.metadata.Clusters)
            if pseudotime:
                
                plot_data['Pseudotime'] = list(self.metadata.Pseudotime)
        
            # Subsetting cells
            plot_data = plot_data.loc[cell_filter,]
        
            # Finding cluster centers for plotting branches
            centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
        
            # Plotting
            if pseudotime:
                
                plt.figure(figsize=(5.5, 5))
                plot_data.sort_values(by="Pseudotime", axis = 0, ascending = True, inplace = True)
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Pseudotime', palette='viridis', marker='.', s=dot_size, linewidth=0)
                norm = plt.Normalize(0, 100)
                sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
                sm.set_array([])
                plt.colorbar(sm)
                plt.legend().remove()
                
                for br in self.branches:
                    
                    cluster_indexes = [clusters.index(cl) for cl in br]
                    x, y = centers[cluster_indexes, 0], centers[cluster_indexes, 1]
                    plt.plot(x, y, markersize=3, marker='o', linewidth=1, color="red", solid_capstyle='round')
                
                plt.xlabel('Umap 1')
                plt.ylabel('Umap 2')
                
                plt.savefig("Trajectories_UMAP-{}_Pseudotime.png".format(n), bbox_inches='tight', dpi=300)

            else:
                
                plt.figure(figsize=(5, 5))
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Clusters', palette=color_palette, marker='.', s=dot_size, linewidth=0)
                legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Clusters', ncol=legend_cols)

                for br in self.branches:
                
                    cluster_indexes = [clusters.index(cl) for cl in br]
                    x, y = centers[cluster_indexes, 0], centers[cluster_indexes, 1]
                    plt.plot(x, y, markersize=3, marker='o', linewidth=1, color="black", solid_capstyle='round')

                plt.xlabel('Umap 1')
                plt.ylabel('Umap 2')
            
                plt.savefig("Trajectories_UMAP-{}_NoPseudotime.png".format(n), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
                
            plt.close()

### ------------------MAIN------------------ ###

import h5py
import igraph
import numpy as np
import pandas as pd
import pickle as pk
import random
import seaborn
import umap

from matplotlib import pyplot as plt
from os import listdir
from scipy.stats import mannwhitneyu, nbinom
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from statsmodels.discrete.discrete_model import NegativeBinomial
from statsmodels.stats.multitest import fdrcorrection
