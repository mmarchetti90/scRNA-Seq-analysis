#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import gc
import gzip
import h5py
import igraph
import numpy as np
import pandas as pd
import pickle as pk
import random
import scanorama
import seaborn
import umap

from datetime import datetime
from itertools import combinations
from matplotlib import pyplot as plt
from os import listdir, makedirs
from os.path import abspath, isdir, exists
from matplotlib import patches as mpatches
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from numpy.typing import NDArray
from scipy.sparse import csr_array, dok_array, load_npz, save_npz, vstack
from scipy.stats import hypergeom, mannwhitneyu, zscore
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, silhouette_score
from sklearn.neighbors import kneighbors_graph
from statsmodels.stats.multitest import fdrcorrection

### CLASSES -------------------------------- ###

class IntegratedAnalysis:

    """
    DESCRIPTION:
    Class to read in, integrate, and analyze published scRNA-Seq datasets (may be raw counts or normalized in different ways)

    NOTES:
    Data must be in one of two formats (note that text files can be gzip compressed):
        - Tab-separated table with genes as rows and cells as columns. First column must be gene symbol.
        - Sparse matrix split into barcodes, features, and matrix files contained in the same folder.
        N.B. The barcodes file first column must be the list of cell_ids/barcodes;
        N.B. The second column of the features file must be gene symbols, but the file can also be a single column of gene symbols;

    TYPICAL PIPELINE:
        0 - load_manifest;
        1 - preprocess_samples;
        2 - find_common_most_variable_features;
        3 - integrate_data;
        4 - merge_integrated_datasets;
        5 - reduce_dimensions;
        6 - cluster_cells;
        7 - find_clusters_markers;
        8 - score_cluster_identity;
        9 - create_mst;
        10 - generate_pseudotime;
        11 - score_cell_cycle or score_gene_set;
    """

    def __init__(self, max_cells_in_memory: int=100000, desired_feature_type: str='Gene Expression', repo_name: str='scrna_analysis') -> None:
            
        # Creating color palette for up to 60 clusters
        """
        # Use the following if there's a lot of clusters
        cluster_colors = [(r/3, g/3, b/3) for r in range(4) for g in range(4) for b in range(4)][1:-1]
        cluster_colors = [cluster_colors[i1 + i2] for i2 in range(4) for i1 in range(0, len(cluster_colors), 4) if i1 + i2 < len(cluster_colors)] # Spacing colors
        random_indexes = [33, 18, 41, 43, 34, 8, 56, 31, 16, 25, 2, 7, 46, 36, 5, 22, 11, 26, 53, 9, 14, 19, 58, 23, 1, 42, 55, 61, 35, 49, 6, 15, 44, 30, 47, 38, 37, 54, 60, 13, 51, 39, 48, 50, 3, 10, 21, 52, 24, 57, 28, 27, 32, 17, 20, 40, 12, 29, 59, 45]
        cluster_colors = [cluster_colors[ri] for ri in random_indexes]
        self.cluster_colors = seaborn.color_palette(cluster_colors, as_cmap=True)
        """
        self.cluster_colors = seaborn.color_palette("tab20") * 10
            
        # Setting maximum number of cells from sparse matrices to load in memory
        self.max_cells_in_memory = max_cells_in_memory
            
        # Setting the selected feature type ('Gene Expression' for RNA data or 'Peak' for ATAC data)
        if desired_feature_type not in ['Gene Expression', 'Peak']:
            self.log_trace("Error: Desired feature type must be 'Gene Expression' or 'Peak'")
            self.log_trace("Defaulting to 'Gene Expression'")
            self.desired_feature_type = 'Gene Expression'
        else:
            self.desired_feature_type = desired_feature_type

        # Define main repository path
        self.repo_name = repo_name

        # Init paths
        self.init_paths()

        # Check if repo exists and load it if it does
        if exists(self.repo_name):
            self.load_repo()
            self.log_trace('Repo was fully loaded')
        else:
            makedirs(self.repo_name)
            self.log_trace('Repo was initialized')

    ### ------------------------------------ ###

    def init_paths(self) -> None:

        # Init repo dirs standard paths
        self.preprocessed_dir = f'{self.repo_name}/preprocessed_data'
        self.hvgs_dir = f'{self.repo_name}/hvgs'
        self.hvgs_file = 'higly_variable_genes.tsv.gz'
        self.integrated_data_dir = f'{self.repo_name}/integrated_data'
        self.integrated_pca_embdeddings_file = 'integrated_pca_embdeddings.tsv.gz'
        self.merged_data_dir = f'{self.repo_name}/merged_data'
        self.lower_dimensions_dir = f'{self.repo_name}/lower_dimensions'
        self.clusters_dir = f'{self.repo_name}/clustering'
        self.cluster_embeddings_file = 'cluster_embeddings.tsv.gz'
        self.cluster_markers_file = 'cluster_markers.tsv.gz'
        self.cluster_identities_file = 'cluster_identities.tsv.gz'
        self.trajectories_dir = f'{self.repo_name}/trajectories'
        self.trajectories_file = 'trajectory_analysis.txt'
        self.pseudotime_file = 'pseudotime.tsv.gz'
        self.cell_cycle_dir = f'{self.repo_name}/cell_cycle'
        self.cell_cycle_scores_file = 'cell_cycle_scores.tsv.gz'
        self.cell_interactions_dir = f'{self.repo_name}/cell_interactions'
        self.cell_interactions_file = 'cell_interactions.tsv.gz'
        self.plot_dir = f'{self.repo_name}/data_plots'

    ### ------------------------------------ ###
    ### TRACE LOG                            ###
    ### ------------------------------------ ###

    def log_trace(self, trace_message: str, timestamp: bool=True) -> None:

        """
        Logs a time-stamped trace message to the console and appends it to the log

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

    ### ------------------------------------ ###
    ### DATA I/O                             ###
    ### ------------------------------------ ###

    def load_manifest(self, manifest_path: str) -> None:
                
        try:
            self.manifest = pd.read_csv(manifest_path, sep='\t')[['data_label', 'dataset_path', 'normalization']]
        except:
            self.log_trace('ERROR: could not read manifest file.')
            self.manifest = pd.DataFrame({'data_label' : [], 'dataset_path' : []})

        if len(np.unique(self.manifest['data_label'].tolist())) < self.manifest.shape[0]:
            self.log_trace('ERROR: labels in data manifest are not unique!')
            self.manifest = pd.DataFrame({'data_label' : [], 'dataset_path' : []})
    
    ### ------------------------------------ ###

    def load_repo(self) -> None:

        """
        Wrapper to load all identifiable elements of a repo
        """

        # Load count merged matrix
        try:
            self.load_merged_datasets(self.merged_data_dir)
        except:
            pass
        
        # Load HVGS
        try:
            self.load_hvgs(f'{self.hvgs_dir}/{self.hvgs_file}')
        except:
            pass

        # Load clusters
        try:
            self.load_clusters(f'{self.clusters_dir}/{self.cluster_embeddings_file}')
        except:
            pass

        # Load cluster markers
        try:
            self.load_cluster_markers(f'{self.clusters_dir}/{self.cluster_markers_file}')
        except:
            pass

        # Load cluster identities
        try:
            self.load_cluster_identities(f'{self.clusters_dir}/{self.cluster_identities_file}')
        except:
            pass

        # Load cell cycle scores
        try:
            self.load_cell_cycle_scores(f'{self.cell_cycle_dir}/{self.cell_cycle_scores_file}')
        except:
            pass
        
        # Load trajectories branches
        try:
            self.load_trajectory_branches(f'{self.trajectories_dir}/{self.trajectories_file}')
        except:
            pass

        # Load pseudotime
        try:
            self.load_trajectory_pseudotime(f'{self.trajectories_dir}/{self.pseudotime_file}')
        except:
            pass
        
        # Load cell_interactions
        try:
            self.load_cell_interactions(f'{self.cell_interactions_dir}/{self.cell_interactions_file}')
        except:
            pass
    
    ### ------------------------------------ ###
    
    def load_all_matrices_in_path(self, data_path: str) -> tuple[list[str], list[csr_array], list[list[str]], list[list[str]]]:
        
        # Get list of matrices to load
        matrices_list = {matrix.split('/')[-1] : f'{data_path}/{matrix}' for matrix in listdir(data_path) if isdir(f'{data_path}/{matrix}')}
        
        # Check that directories contain the needed files
        needed_files = ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']
        matrices_list = {key : value for (key, value) in matrices_list.items() if sum([nf in file for nf in needed_files for file in listdir(value)]) >= 3}
        
        # Get datasets labels
        datasets_labels = list(matrices_list.keys())
        
        # Load datasets as sparse matrices and genes lists
        datasets, cells_list, genes_list = [], [], []
        for (label, path) in matrices_list.items():
            # Load matrix as sparse file of cells-by-genes
            barcodes, features, sparse_matrix = self.load_sparse_matrix(path, save_npz_if_missing=False)
            # Store data
            datasets.append(sparse_matrix)
            cells_list.append(barcodes)
            genes_list.append(features)

        return datasets_labels, datasets, cells_list, genes_list

    ### ------------------------------------ ###
    
    def load_cell_cycle_scores(self, data_path) -> None:
        
        try:
            self.cell_cycle_scores = pd.read_csv(data_path, sep='\t', header=0)
        except:
            self.log_trace("ERROR: couldn't read cell cycle file")
            self.cell_cycle_scores = pd.DataFrame()

    ### ------------------------------------ ###

    def load_cell_interactions(self, data_path) -> None:

        try:
            self.cell_interactions = pd.read_csv(data_path, sep='\t', header=0)
        except:
            self.log_trace("ERROR: couldn't read cell interactions file")
            self.cell_interactions = pd.DataFrame()
    
    ### ------------------------------------ ###
    
    def load_clusters(self, data_path: str) -> None:
        
        try:
            self.clusters = pd.read_csv(data_path, sep='\t', header=None).to_numpy().flatten()
            self.clusters_dir = '/'.join(data_path.split('/')[:-1])
        except:
            self.log_trace("ERROR: couldn't read clusters annotation file")
            self.clusters = np.array([])
            self.clusters_dir = ''
    
    ### ------------------------------------ ###
    
    def load_cluster_identities(self, data_path: str) -> None:
        
        try:
            self.cluster_identities = pd.read_csv(data_path, sep='\t', header=0)
        except:
            self.log_trace("ERROR: couldn't read clusters identity file")
            self.cluster_identities = pd.DataFrame()
    
    ### ------------------------------------ ###
    
    def load_cluster_markers(self, data_path: str) -> None:
        
        try:
            self.cluster_markers = pd.read_csv(data_path, sep='\t', header=0)
        except:
            self.log_trace("ERROR: couldn't read clusters markers file")
            self.cluster_markers = pd.DataFrame()
    
    ### ------------------------------------ ###
        
    def load_hvgs(self, data_path: str) -> None:
            
        try:
            self.hvgs = pd.read_csv(data_path, sep='\t', header=None).to_numpy().flatten().tolist()
        except:
            self.log_trace("ERROR: couldn't read HVGs file")
            self.hvgs = []
    
    ### ------------------------------------ ###
    
    def load_merged_datasets(self, data_path: str) -> None:
        
        try:
            self.datasets_metadata = pd.read_csv(f'{data_path}/datasets_metadata.tsv.gz', sep='\t', header=0)
            self.all_cells, self.all_genes, self.all_data = self.load_sparse_matrix(data_path, save_npz_if_missing=False)
            self.merged_data_dir = data_path
        except:
            self.log_trace("ERROR: couldn't read merged dataset files")
            self.merged_data_dir = ''
    
    ### ------------------------------------ ###
    
    def load_reduced_dimensions(self, lower_dimensions_dir) -> tuple[PCA, np.ndarray, int, umap.UMAP, np.ndarray]:
        
        try:
            # Find necessary items in lower_dimensions_dir
            pca_model_path = [f'{lower_dimensions_dir}/{file}' for file in listdir(lower_dimensions_dir) if 'pca_optimalcomp-' in file and file.endswith('.pkl')]
            pca_data_path = [f'{lower_dimensions_dir}/{file}' for file in listdir(lower_dimensions_dir) if file == 'pca_data.tsv.gz'][0]
            umap_model_path = [f'{lower_dimensions_dir}/{file}' for file in listdir(lower_dimensions_dir) if 'umap_neighbors-' in file and file.endswith('.pkl')][0]
            umap_data_path = [f'{lower_dimensions_dir}/{file}' for file in listdir(lower_dimensions_dir) if file == 'umap_data.tsv.gz'][0]
            
            # Load PCA model
            if len(pca_model_path) > 0:
                pca_model = pk.load(open(pca_model_path[0], "rb"))
                # Extract PCA optimal components
                x0, x1 = 0, len(pca_model.explained_variance_) - 1
                y0, y1 = pca_model.explained_variance_[x0], pca_model.explained_variance_[x1]
                gradient = (y1 - y0) / (x1 - x0)
                intercept = y0
                difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(pca_model.explained_variance_)]
                optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
            else:
                pca_model = PCA()
                optimal_pca_components = 1000000
            
            # Load UMAP model
            umap_model = pk.load(open(umap_model_path, "rb"))
            
            # Load pca_data and umap_data
            pca_data = pd.read_csv(pca_data_path, sep='\t', header=None).to_numpy()[:, :optimal_pca_components]
            umap_data = pd.read_csv(umap_data_path, sep='\t', header=None).to_numpy()
            
            return pca_model, pca_data, optimal_pca_components, umap_model, umap_data
        
        except:
            
            self.log_trace('ERROR: make sure all the needed files exist in folder.')
            return PCA(), np.array([]), 0, umap.UMAP(), np.array([])
    
    ### ------------------------------------ ###
    
    def load_sparse_matrix(self, main_dir: str, save_npz_if_missing: bool=True) -> tuple[list[str], list[str], csr_array]:
        
        try:
            # Find necessary files
            barcodes_path, features_path, matrix_path = [l for l in listdir(main_dir) if 'barcodes' in l][0], [l for l in listdir(main_dir) if 'features' in l][0], [l for l in listdir(main_dir) if 'matrix' in l]
            
            # Read cell_identities file
            barcodes = self.read_barcodes_file(f'{main_dir}/{barcodes_path}')
            
            # Read genes list file
            features, features_type = self.read_features_file(f'{main_dir}/{features_path}')
            
            # Get number of barcodes and features
            cells_num, genes_num = len(barcodes), len(features)
            
            # Check if npz sparse matrix exists
            npz_file = [l for l in listdir(main_dir) if l.endswith('.npz')]
            if len(npz_file):
                sparse_matrix = load_npz(f'{main_dir}/{npz_file[0]}')
                # Filter features based on type
                if len(features_type):
                    features_filter = np.where(np.array(features_type) == self.desired_feature_type)[0]
                    features = np.array(features)[features_filter].tolist()
                    sparse_matrix = sparse_matrix[:, features_filter]
            else:
                # Constructing a new csr_array
                matrix_path = matrix_path[0]
                
                # Check matrix files, finding best way to break it into manageable chunks of self.max_cells_in_memory
                matrix_separator, matrix_breaks = self.check_matrix_file(f'{main_dir}/{matrix_path}', cells_num, genes_num, self.max_cells_in_memory)
            
                # Load matrix as a cells-by-genes scipy.sparse.csr_array
                sparse_matrix = dok_array((cells_num, genes_num))
                for start, stop in zip(matrix_breaks[:-1], matrix_breaks[1:]):
                    # Load matrix chunk
                    matrix_chunk = self.read_matrix_chunk(f'{main_dir}/{matrix_path}', start, stop, matrix_separator)
                    # Split into data (normalized counts), rows (cells), columns (genes)
                    data = matrix_chunk[:, 2]
                    rows = matrix_chunk[:, 1].astype(int) - 1
                    columns = matrix_chunk[:, 0].astype(int) - 1
                    # Update matrix
                    sparse_matrix[rows, columns] = data

                # Filter features based on type
                if len(features_type):
                    features_filter = (np.array(features_type) == self.desired_feature_type)
                    features = np.array(features)[features_filter].tolist()
                    sparse_matrix = sparse_matrix[:, features_filter]

                # Convert to csr_array
                sparse_matrix = csr_array(sparse_matrix)
                
                # Save as npz file
                if save_npz_if_missing:
                    save_npz(f'{main_dir}/{matrix_path}.npz', sparse_matrix, compressed=True)
            return barcodes, features, sparse_matrix
        except:
            return [], [], csr_array((0, 0))

    ### ------------------------------------ ###
    
    def load_trajectory_branches(self, data_path: str) -> None:
        
        try:
            self.trajectories_dir = '/'.join(abspath(data_path).split('/')[:-1])
            self.trajectory_file_name = data_path.split('/')[-1]
            self.branches = [[int(cluster) for cluster in line.split('\t')] for line in open(data_path).read().split('\n')]
        except:
            self.log_trace("ERROR: couldn't read trajectory branches file")
            self.trajectories_dir = '/'.join(abspath(data_path).split('/')[:-1])
            self.trajectory_file_name = ''
            self.branches = []
    
    ### ------------------------------------ ###
    
    def load_trajectory_pseudotime(self, data_path: str) -> None:
        
        try:
            self.pseudotimes = pd.read_csv(data_path, sep='\t', header=None).iloc[:, 0].to_numpy()
        except:
            self.log_trace("ERROR: couldn't read trajectory pseudotime file")
            self.pseudotimes = np.array([])

    ### ------------------------------------ ###
        
    def merge_datasets(self, datasets_labels: list[str], datasets: list[csr_array], cells_list: list[list[str]], genes_list: list[list[str]], desired_genes: list[str]=[]) -> tuple[pd.DataFrame, list[str], list[str], csr_array]:
        
        # Create a metadata dataframe to store dataset labels and number of cells
        metadata = pd.DataFrame({
            'datasets_labels' : datasets_labels,
            'cell_num' : [len(cl) for cl in cells_list]
        })
        
        # Merge barcodes
        all_barcodes = np.concatenate(cells_list)
        all_barcodes = all_barcodes.tolist()
        
        # Use specified list of genes or find common genes if not specified
        if len(desired_genes) > 0:
            # Make sure gene list is valid
            all_genes = np.array(desired_genes)
            for gl in genes_list:
                all_genes = all_genes[np.isin(all_genes, gl)]
            if not len(all_genes):
                self.log_trace('ERROR: list not valid')
                return metadata, [], [], csr_array((0,0))
        else:
            # Use genes in common among datasets
            all_genes = np.array([])
            for _,gl in enumerate(genes_list):
                if not len(all_genes):
                    all_genes = np.array(gl.copy())
                else:
                    all_genes = all_genes[np.isin(all_genes, gl)]
        all_genes = all_genes.tolist()
        
        # Merging datasets
        all_data = []
        for (ds, gl) in zip(datasets, genes_list):
            # Find indexes of common_genes
            genes_index = [gl.index(g) for g in all_genes]
            # Store matrix subset
            all_data.append(csr_array(ds[:, genes_index]))
        if len(all_data):
            all_data = csr_array(vstack(all_data))
        else:
            all_data = csr_array((0, 0))
        
        return metadata, all_barcodes, all_genes, all_data
    
    ### ------------------------------------ ###

    def save_sparse_matrix_to_path(self, main_path: str, label: str, barcodes: list[str], features: list[str], matrix: csr_array, save_mtx: bool=True) -> str:
        
        # Saves a sparse matrix in nzip and 10X Genomics format to a desired path after preprocessing
        # N.B. 10X Genomics mtx matrix is saved only if save_mtx = True
        # N.B. The 10X Genomics format was chosen since it's quite human readable and can easily be
        # imported in R Seurat if one desires
        
        # Create dir to store sparse matrix
        if not exists(f'{main_path}/{label}'):
            makedirs(f'{main_path}/{label}')
        
        # Save barcodes
        out_text = ('\n'.join(barcodes)).encode()
        out_text += b'\n'
        with gzip.open(f'{main_path}/{label}/{label}_barcodes.tsv.gz', 'wb') as output:
            output.write(out_text)
        
        # Save features
        out_text = ('\n'.join(features)).encode()
        out_text += b'\n'
        with gzip.open(f'{main_path}/{label}/{label}_features.tsv.gz', 'wb') as output:
            output.write(out_text)
        
        # Save matrix as npz
        save_npz(f'{main_path}/{label}/{label}_matrix.mtx.npz', matrix, compressed=True)

        # Save matrix as 10X Genomics mtx file
        if save_mtx:
            self.csr_to_mtx(matrix, out_name=f'{main_path}/{label}/{label}_matrix.mtx.gz', max_cells=self.max_cells_in_memory)
        
        return f'{main_path}/{label}'

    ### ------------------------------------ ###
    
    @staticmethod
    def check_matrix_file(path: str, cells_num: int, genes_num: int, max_cells: int) -> tuple[str, list[int]]:
        
        # Read matrix to find how many lines to skip when loading as pandas as well as the optimal
        # way to break the matrix into chunks of max_cells_in_memory cells
        # N.B. Only checks the first N_lines_to_check lines for headers
        
        # Open matrix file
        if path.endswith('.gz'):
            matrix_file = gzip.open(path,'rb')
        else:
            matrix_file = open(path, 'r')
        
        # Parse matrix file
        N_lines_to_check = 10
        header_lines, breaks = 0, []
        previous_line_cell_index, cell_count = 0, 0
        for (index, line) in enumerate(matrix_file):
            if not len(line): # Skip empty line at bottom of text, if any
                continue
            try:
                if isinstance(line, bytes):
                    line = line.decode('utf8')
                separator = '\t' if str('\t') in line else ' '
                r, c, n = line.split(separator)
                r, c, n = int(r), int(c), float(n) # Columns are 1-N so no need to c + 1 to account for column with gene names
                if (r == genes_num or c == cells_num) and index < N_lines_to_check:
                    header_lines += 1
                else:
                    if c != previous_line_cell_index:
                        previous_line_cell_index = c
                        cell_count += 1
                if cell_count == max_cells:
                    breaks.append(index + 1)
                    cell_count = 0
            except:
                header_lines += 1
        breaks.append(index + 2)
        
        # Adding header lines to breaks
        breaks = [header_lines] + breaks
        
        # Closing matrix file
        matrix_file.close()
        
        return separator, breaks
    
    ### ------------------------------------ ###
    
    @staticmethod
    def read_barcodes_file(path: str) -> list[str]:
        
        barcodes = pd.read_csv(path, header=None,  sep='\t').iloc[:, 0].to_list()
        if barcodes[0].lower() in ['barcode', 'barcodes', 'cell', 'cell_id']:
            barcodes = barcodes[1:]
        
        return barcodes
    
    ### ------------------------------------ ###
    
    @staticmethod
    def read_features_file(path: str) -> tuple[list[str], list[str]]:
        
        features = pd.read_csv(path, header=None,  sep='\t')
        if features.shape[1] == 1:
            features_type = []
            features = features.iloc[:, 0].to_list()
        else:
            features_type = features.iloc[:, 2].to_list()
            features = features.iloc[:, 1].to_list()
        if features[0].lower() in ['gene', 'symbol', 'gene_symbol', 'genesymbol', 'gene_name', 'genename']:
            features = features[1:]
            if len(features_type):
                features_type = features_type[1:]
        
        return features, features_type
    
    ### ------------------------------------ ###
    
    @staticmethod
    def read_matrix_chunk(matrix_path: str, start: int, stop: int, matrix_separator: str) -> np.ndarray:
        
        # Reads a chunk of a 10X Genomics formatted mtx file
        
        # Define rows to skip and how many to read
        skip_rows = start
        n_rows = stop - start
        
        # Load count matrix
        matrix_chunk = pd.read_csv(matrix_path, header=None, skiprows=skip_rows, nrows=n_rows, skip_blank_lines=True, sep=matrix_separator).to_numpy()
        
        return matrix_chunk

    ### ------------------------------------ ###
    ### PREPROCESSING                        ###
    ### ------------------------------------ ###
    
    def preprocess_samples(self, min_cell_raw_counts: int=500, min_detected_genes: int=500):
        
        self.log_trace('Preprocessing data')
        
        # Create dir to store sparse matrices after preprocessing
        if not exists(self.preprocessed_dir):
            makedirs(self.preprocessed_dir)
        
        # Load and normalize samples, then store them in a sparse matrix
        for _,row in self.manifest.iterrows():
            # Extract info
            label, path, norm_type = row
            self.log_trace(f'  * Reading {path}', False)
            
            # Load data
            if isdir(path): # Data is a sparse matrix
                barcodes, features, matrix = self.load_sparse_matrix(path, save_npz_if_missing=False)
            elif path.endswith('.h5'): # Data in h5 format
                barcodes, features, matrix = self.h5_to_csr(path, self.desired_feature_type)
            else: # Data is a full matrix
                data = pd.read_csv(path, sep='\t', skip_blank_lines=True)
                self.log_trace(f'    * Converting {path} to sparse matrix', False)
                barcodes, features, matrix = self.full_to_sparse(data)
                del data
                gc.collect()
            
            # Skip if dataframe could not be read
            if len(barcodes) == 0:
                self.log_trace(f'    * WARNING: could not read file {path}', False)
                continue
            
            # Remove cells with number of genes detected < min_detected_genes
            self.log_trace('    * Removing cells with too few genes detected', False)
            barcodes, matrix = self.remove_cells_with_few_genes(barcodes, matrix, min_detected_genes)
            
            # Normalize data
            if norm_type == 'none': # Performing PFlog1pPF normalization
                # Remove cells with read count < min_cell_raw_counts
                self.log_trace('    * Removing cells with too few counts', False)
                barcodes, matrix = self.remove_cells_with_low_counts(barcodes, matrix, min_cell_raw_counts)
                # Normalize data
                self.log_trace('    * Normalizing counts', False)
                matrix = self.normalize_data(matrix, skip_first_proportional_fitting=False)
            elif norm_type in ['cpm', 'rpkm']: # Performing log1p normalization and proportional fitting
                self.log_trace('    * Normalizing counts', False)
                matrix = self.normalize_data(matrix, skip_first_proportional_fitting=True)
            elif norm_type in ['pflog1ppf', 'pool_norm']: # Preferred normalizations, data is passed as is
                pass
            else:
                self.log_trace(f'    * WARNING: could not determine normalization method of {path}', False)
                self.log_trace('    * Data will be used as is...', False)
            
            # Append label to cell names
            barcodes = [b + '_' + label for b in barcodes]
            
            # Saving sparse matrix to temporary file
            self.log_trace('    * Saving sparse matrix', False)
            _ = self.save_sparse_matrix_to_path(self.preprocessed_dir, label, barcodes, features, matrix)
        
        self.log_trace('Preprocessing done!')
        self.log_trace(f'  * Data was saved to sparse matrices in {self.preprocessed_dir}', False)

    ### ------------------------------------ ###
        
    def normalize_data(self, mtx: csr_array, skip_first_proportional_fitting: bool=False) -> csr_array:
        
        # Normalizes gene counts by proportional fitting + log1p + proportional fitting
        # N.B. first proportional fitting should be skipped if data was already CPM or RPKM normalized
        
        # First proportional fitting
        if not skip_first_proportional_fitting:
            mtx = self.proportional_fitting(mtx)
        
        # Log1p transformation
        mtx.data = np.log1p(mtx.data)
        
        # Second proportional fitting
        mtx = self.proportional_fitting(mtx)
        
        return mtx
    
    ### ------------------------------------ ###
    
    @staticmethod
    def csr_to_mtx(csr: csr_array, out_name: str, max_cells: int=1000) -> None:
        
        # Save csr_array as 10X Genomics mtx file (100 cells at a time to save memory)
        with gzip.open(out_name, 'wb') as mtx_out:
            # Write header
            out_text = '\n'.join([
                "%%MatrixMarket matrix coordinate numeric general",
                f'{csr.shape[1]}\t{csr.shape[0]}\t{csr.count_nonzero()}'
            ])
            out_text = out_text.encode()
            out_text += b'\n'
            mtx_out.write(out_text)
            for i in range(csr.shape[0] // max_cells + 1):
                # Creating indexes for mtx rows (genes) and columns (cells)
                # Note that in the csr, rows and columns are transposed
                cols, rows = np.where(csr_array(csr[i * max_cells : (i + 1) * max_cells,:]).toarray() != 0)
                cols += i * max_cells
                # Get values
                counts = np.array(csr[cols, rows]).flatten()
                # Adjust columns and rows indices (mtx files are 1-based)
                rows += 1
                cols += 1
                # Format text
                out_text = '\n'.join([f'{r}\t{c}\t{n}' for r, c, n in zip(rows, cols, counts)])
                out_text = out_text.encode()
                out_text += b'\n'
                # Write text
                mtx_out.write(out_text)
    
    ### ------------------------------------ ###
    
    @staticmethod
    def full_to_sparse(full: pd.DataFrame) -> tuple[list[str], list[str], csr_array]:
        
        # Converts a full matrix to a csr_array
        
        # Get features and barcodes
        barcodes = full.columns[1:].to_list()
        features = full['GeneSymbol'].to_list()
        
        # Create counts matrix
        matrix = csr_array(full.iloc[:, 1:].T.to_numpy())
        
        return barcodes, features, matrix
    
    ### ------------------------------------ ###
    
    @staticmethod
    def h5_to_csr(h5_path: str, desired_feature_type: str='') -> tuple[list[str], list[str], csr_array]:

        # Reading HDF5 data
        with h5py.File(h5_path, 'r') as h5_data:
            # Unique cell identifiers
            barcodes = np.array(h5_data['matrix/barcodes'], dtype=str).tolist()
            cells_num = len(barcodes)
                    
            # Genes
            #gene_ids = np.array(h5_data['matrix']['features']['id']).astype(str)
            gene_names = np.array(h5_data['matrix/features/name'], dtype=str)
            genes_num = len(gene_names)
            try:
                features_type = np.array(h5_data['matrix/features/feature_type']).astype(str)
            except:
                features_type = []
            
            # Init sparse matrix
            sparse_matrix = dok_array((cells_num, genes_num))
            
            # Parse sparse
            for N,(start,stop) in enumerate(zip(np.array(h5_data['matrix/indptr'], dtype=int), np.array(h5_data['matrix/indptr'], dtype=int)[1:])):
                # Extracting cell data
                counts = np.array(h5_data['matrix/data'], dtype=int)[start : stop]
                indexes = np.array(h5_data['matrix/indices'], dtype=int)[start : stop]
                # Update matrix
                sparse_matrix[N, indexes] = counts
            
            # Filter features based on type
            if len(features_type) and desired_feature_type in ['Gene Expression', 'Peak']:
                features_filter = np.where(np.array(features_type) == desired_feature_type)[0]
                features = np.array(gene_names)[features_filter].tolist()
                sparse_matrix = sparse_matrix[:, features_filter]
    
            # Convert to csr_array
            sparse_matrix = csr_array(sparse_matrix)
        
        return barcodes, features, sparse_matrix
    
    ### ------------------------------------ ###
        
    @staticmethod
    def proportional_fitting(mtx: csr_array, target_size: float=-1) -> csr_array:
        
        # Get total cell counts
        total_counts = np.asarray(mtx.sum(axis=1)).reshape(-1)
        
        # Compute correction factors
        if target_size == -1:
            target_size = np.median(total_counts).astype(float)
        else:
            pass
        correction_factors = target_size / total_counts
        correction_factors = correction_factors.reshape((len(correction_factors), 1))
        
        # Correct matrix
        mtx = csr_array(mtx.multiply(correction_factors))

        return mtx
    
    ### ------------------------------------ ###
    
    @staticmethod
    def remove_cells_with_few_genes(brcs, mtx, min_genes):
        
        # Find cells with less than min_cnts total raw counts from the count matrix
        good_cells = np.asarray((mtx != 0).sum(axis=1) >= min_genes).reshape(-1)

        # Removes bad cells from barcodes list and matrix
        brcs = np.array(brcs)[good_cells].tolist()
        mtx = mtx[good_cells, ]
        
        return brcs, mtx
    
    ### ------------------------------------ ###
    
    @staticmethod
    def remove_cells_with_low_counts(brcs, mtx, min_cnts):
        
        # Find cells with less than min_cnts total raw counts from the count matrix
        good_cells = np.asarray(mtx.sum(axis=1) >= min_cnts).reshape(-1)

        # Removes bad cells from barcodes list and matrix
        brcs = np.array(brcs)[good_cells].tolist()
        mtx = mtx[good_cells, ]
        
        return brcs, mtx

    ### ------------------------------------ ###
    ### HIGHLY VARIABLE GENES                ###
    ### ------------------------------------ ###
    
    def find_common_highly_variable_genes(self, data_path: str, x_low_cutoff: float=0.1, x_high_cutof: float=8, y_low_cutoff: float=1, y_high_cutoff: float=1e9, nbins: int=20, max_features: int=3000, union_mode: bool=False):
        
        # Wrapper for finding highly variable genes (HVGs) common to all datasets
        # If union_mode = True, all HVGs from individual datasets are combined
        self.log_trace(f'Finding {"all" if union_mode else "common"} highly variable genes')
        
        # Create dir to store HVGs
        if not exists(self.hvgs_dir):
            makedirs(self.hvgs_dir)
        
        # Load matrices
        datasets_labels, datasets, _, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Find variable features in each dataset, then intersect
        hvgs = np.array([])
        for label, ds, gl in zip(datasets_labels, datasets, genes_list):
            self.log_trace(f'  * Finding highly variable genes for {label}', False)
            # Finding HVGs indexes
            dataset_hvgs = self.find_highly_variable_genes(ds, x_low_cutoff, x_high_cutof, y_low_cutoff, y_high_cutoff, nbins, max_features)
            # Extracting gene names
            gene_names = np.array(gl)[dataset_hvgs]
            # Storing data
            if not len(hvgs):
                hvgs = gene_names.copy()
            else:
                if not union_mode:
                    hvgs = hvgs[np.isin(hvgs, gene_names)]
                else:
                    hvgs = np.unique(np.concatenate([hvgs, gene_names]))
            self.log_trace(f'    * Found {len(dataset_hvgs)} highly variable genes', False)
        
        # Store HVGs
        self.hvgs = list(hvgs)
        
        # Save HVGs to file
        out_text = ('\n'.join(self.hvgs)).encode()
        out_text += b'\n'
        with gzip.open(f'{self.hvgs_dir}/{self.hvgs_file}', 'wb') as output:
            output.write(out_text)
        self.log_trace(f'Found {len(hvgs)} common highly variable genes')
        self.log_trace(f'  * Highly variable genes were saved to {self.hvgs_dir}/{self.hvgs_file}', False)
    
    ### ------------------------------------ ###
    
    def find_highly_variable_genes_in_whole_dataset(self, data_path: str, x_low_cutoff: float=0.1, x_high_cutof: float=8, y_low_cutoff: float=1, y_high_cutoff: float=1e9, nbins: int=20, max_features: int=3000):
        
        # Find highly variable genes after merging all datasets into one
        # Use this if find_common_highly_variable_genes returns very few genes
        self.log_trace('Finding highly variable genes in merged datasets')
        
        # Create dir to store HVGs
        if not exists(self.hvgs_dir):
            makedirs(self.hvgs_dir)
        
        # Load matrices
        _, datasets, _, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Merge data
        _, _, all_genes, all_data = self.merge_datasets([], datasets, [[], []], genes_list, desired_genes=[])
        
        # Finding HVGs indexes
        hvgs = self.find_highly_variable_genes(all_data, x_low_cutoff, x_high_cutof, y_low_cutoff, y_high_cutoff, nbins, max_features)
        
        # Extracting gene names
        hvgs = np.array(all_genes)[hvgs]
    
        # Store HVGs
        self.hvgs = list(hvgs)
        
        # Save HVGs to file
        out_text = ('\n'.join(self.hvgs)).encode()
        out_text += b'\n'
        with gzip.open(f'{self.hvgs_dir}/{self.hvgs_file}', 'wb') as output:
            output.write(out_text)
        self.log_trace(f'Found {len(hvgs)} highly variable genes')
        self.log_trace(f'  * Highly variable genes were saved to {self.hvgs_dir}/{self.hvgs_file}', False)
    
    ### ------------------------------------ ###
    
    def find_highly_variable_genes(self, mtx: csr_array, x_low_cutoff: float=0.1, x_high_cutof: float=8, y_low_cutoff: float=1, y_high_cutoff: float=1e9, nbins: int=20, max_features: int=3000) -> list[str]:
        
        # Graph based, similar to Seurat's FindVariableGenes
        
        # Transposing matrix from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_array slicing is efficient along rows
        mtx_t = mtx.transpose()
        
        # Calculate features mean and variance
        # Remember that variance = mean of squared values - square of mean values
        self.log_trace('  * Computing genes mean and variance', False)
        mtx_squared = mtx_t.copy()
        mtx_squared.data **= 2
        gene_mean = np.asarray(mtx_t.mean(axis=1)).reshape(-1)
        gene_variance = np.asarray(mtx_squared.mean(axis=1)).reshape(-1) - (gene_mean ** 2)
        
        # Calculating e^mean and log(variance/mean) (X and Y respectively) for each feature
        var_plot = pd.DataFrame({
            'GeneID' : range(0, mtx_t.shape[0]),
            'X' : np.exp(gene_mean),
            'Y' : np.log(gene_variance / (gene_mean + 0.00001))
        })
        min_x, max_x = max(x_low_cutoff, var_plot.X.min()), min(x_high_cutof, var_plot.X.max())
        bin_width = (max_x - min_x) / nbins

        # Binning and computing log(variance/mean) z-scores
        self.log_trace('  * Binning genes and computing log(variance/mean) z-scores', False)
        most_variable_features = {'GeneID' : [], 'ZScores' : []}
        for x_start in np.arange(min_x, max_x, bin_width):
            x_stop = x_start + bin_width
            binned = var_plot.loc[(var_plot['X'] > x_start) & (var_plot['X'] <= x_stop)]
            var_z_scores = (binned.Y - binned.Y.mean()) / binned.Y.std()
            most_variable_features['GeneID'].extend(list(binned.GeneID))
            most_variable_features['ZScores'].extend(list(var_z_scores))
        most_variable_features = pd.DataFrame(most_variable_features)
        most_variable_features.sort_values(by = 'ZScores', axis = 0, ascending = False, inplace = True)

        # Removing features with too low/high z_score, then selecting top ones
        most_variable_features = most_variable_features.loc[
            (most_variable_features.ZScores >= y_low_cutoff) &
            (most_variable_features.ZScores <= y_high_cutoff)
        ].iloc[:max_features,].GeneID.tolist()
        
        return most_variable_features
    
    ### ------------------------------------ ###
    ### DATA INTEGRATION                     ###
    ### ------------------------------------ ###
    
    def integrate_data(self, data_path: str='', use_highly_variable_genes: bool=True, union: bool=False, adjust_library_size: bool=True, use_scanorama_embeddings: bool=True, max_pca_components: int=50) -> None:
        
        # Data integration with Scanorama
        self.log_trace('Integrating data using Scanorama')
        
        # Define path to preprocessed sparce matrices folder
        if not len(data_path):
            try:
                data_path = self.preprocessed_dir
            except:
                self.log_trace('  * ERROR: missing data_path value', False)
                return None
        
        # Check that data_path exists
        if not exists(data_path):
            self.log_trace('  * ERROR: bad data_path value', False)
            return None
        
        # Create directory to store the integrated sparse matrix
        if not exists(self.integrated_data_dir):
            makedirs(self.integrated_data_dir)
        
        # Load matrices for Scanorama
        datasets_labels, datasets, cells_list, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Correcting for library size
        if adjust_library_size:
            datasets_target_size = np.median(np.concatenate([ds.sum(axis=1).ravel() for ds in datasets])).astype(float)
            datasets = [self.proportional_fitting(ds, datasets_target_size) for ds in datasets]
        
        # Filter for common HVGs
        if use_highly_variable_genes:
            # Merging datasets using HVGs, then splitting into individual datasets. This so as to filter datasets for HVGs
            datasets_metadata, _, all_genes, all_data = self.merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=self.hvgs)
            # Splitting all_data into its constituent datasets
            for index in range(datasets_metadata.shape[0]):
                # Get dataset info
                label = datasets_metadata.iloc[index, 0]
                start_row = datasets_metadata.iloc[:index, 1].toarray.sum()
                end_row = datasets_metadata.iloc[:index + 1, 1].toarray.sum()
                datasets[index] = csr_array(all_data[start_row : end_row])
            genes_list = [all_genes for _ in range(len(datasets))]
            
        # Integrate data
        if not use_scanorama_embeddings:
            corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=False, ds_names=datasets_labels, approx=True, knn=20, hvg=None, return_dense=False, union=union, seed=42, batch_size=self.max_cells_in_memory) # type: ignore
        else:
            integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True, ds_names=datasets_labels, approx=True, knn=20, hvg=None, return_dense=False, union=union, seed=42, batch_size=self.max_cells_in_memory) # type: ignore
            integrated = np.concatenate(integrated, axis=0)[:, :max_pca_components]
            self.pca_data = integrated
            pd.DataFrame(integrated).to_csv(f'{self.integrated_data_dir}/{self.integrated_pca_embdeddings_file}', sep='\t', header=False, index=False)
        
        # Save matrices
        for (label, matrix, barcodes) in zip(datasets_labels, corrected, cells_list):
            self.log_trace(f'Saving {label} corrected matrix to file')
            self.save_sparse_matrix_to_path(self.integrated_data_dir, label, barcodes, genes, matrix, save_mtx=True) # type: ignore
        self.log_trace('All done!')
        self.log_trace('  * Data was saved to sparse matrices in {self.integrated_data_dir}', False)
        
        # Merging integrated dataset
        self.merge_integrated_datasets(self.integrated_data_dir, False)
    
    ### ------------------------------------ ###
    
    def merge_integrated_datasets(self, data_path: str='', adjust_library_size: bool=False) -> None:
        
        self.log_trace('Merging integrated data to single dataset')
        
        # Define path to integrated sparce matrices folder
        if not len(data_path):
            try:
                data_path = self.integrated_data_dir
            except:
                self.log_trace('  * ERROR: missing data_path value', False)
                return None
        
        # Check that data_path exists
        if not exists(data_path):
            self.log_trace('  * ERROR: bad data_path value', False)
            return None
        
        # Create directory to store the merged sparse matrix
        if not exists(self.merged_data_dir):
            makedirs(self.merged_data_dir)
        
        # Load integrated data
        datasets_labels, datasets, cells_list, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Correcting for library size
        if adjust_library_size:
            datasets_target_size = np.median(np.concatenate([ds.sum(axis=1).ravel() for ds in datasets])).astype(float) # type: ignore
            datasets = [self.proportional_fitting(ds, datasets_target_size) for ds in datasets]
        
        # Merge datasets
        datasets_metadata, all_cells, all_genes, all_data = self.merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=[])
        
        # Store in class
        self.datasets_metadata = datasets_metadata
        self.all_cells = all_cells
        self.all_genes = all_genes
        self.all_data = all_data
        
        # Save metadata
        self.log_trace('  * Saving datasets metadata to file', False)
        datasets_metadata.to_csv(f'{self.merged_data_dir}/datasets_metadata.tsv.gz', sep='\t', index=False)

        # Save barcodes
        self.log_trace('  * Saving merged barcodes to file', False)
        out_text = ('\n'.join(all_cells)).encode()
        out_text += b'\n'
        with gzip.open(f'{self.merged_data_dir}/all_barcodes.tsv.gz', 'wb') as output:
            output.write(out_text)
        
        # Save features
        self.log_trace('  * Saving merged gene names to file', False)
        out_text = ('\n'.join(all_genes)).encode()
        out_text += b'\n'
        with gzip.open(f'{self.merged_data_dir}/all_features.tsv.gz', 'wb') as output:
            output.write(out_text)
        
        # Save matrix as npz
        self.log_trace('  * Saving integrated sparse matrix to file', False)
        save_npz(f'{self.merged_data_dir}/all_data_matrix.mtx.npz', all_data, compressed=True)

        # Save matrix as 10X Genomics mtx file
        self.log_trace('  * Saving integrated sparse matrix to file in 10X Genomics format', False)
        self.csr_to_mtx(all_data, out_name=f'{self.merged_data_dir}/all_data_matrix.mtx.gz', max_cells=self.max_cells_in_memory)
        
        self.log_trace('All done!')
        self.log_trace('  * Data was saved to sparse matrices in {self.merged_data_dir}', False)

    ### ------------------------------------ ###
    ### DIMENSIONALITY REDUCTION             ###
    ### ------------------------------------ ###
    
    def reduce_dimensions(self, use_integrated_pca: bool=True, integrated_pca_path: str='', use_highly_variable_genes: bool=True, max_pca_components: int=50, neighbors: int=30, n_permutations: int=3) -> None:
        
        # Wrapper for umap_embedding and pca_umap_embdedding functions
        # N.B. PCA/UMAP models and embeddings are not kept in memory, but have to be loaded with self.load_reduced_dimensions afterwards
        if use_integrated_pca:
            self.umap_embedding(integrated_pca_path, neighbors, n_permutations)
        else:
            self.pca_umap_embdedding(use_highly_variable_genes, max_pca_components, neighbors, n_permutations)
    
    ### ------------------------------------ ###
    
    def pca_umap_embdedding(self, use_highly_variable_genes: bool, max_pca_components: int=50, neighbors: int=30, n_permutations: int=3):
        
        self.log_trace('Reducing dataset dimensions')
        
        # Create dir to store dimensionality reduction data
        if not exists(self.lower_dimensions_dir):
            makedirs(self.lower_dimensions_dir)
        
        # Filter for HVGs
        all_data = csr_array(self.all_data.copy())
        if use_highly_variable_genes:
            all_data = csr_array(all_data[:, [self.all_genes.index(hvg) for hvg in self.hvgs]])
            
        # Create subsets of self.max_cells_in_memory cells to fit the PCA and UMAP models to
        # N.B. All cells are passed through 3 times in different groupings
        if all_data.shape[0] < self.max_cells_in_memory:
            cell_subsets = self.make_cell_subsets(all_data.shape[0], self.max_cells_in_memory)
        else:
            cell_subsets = [subset for _ in range(n_permutations) for subset in self.make_cell_subsets(all_data.shape[0], self.max_cells_in_memory)]
            
        # Calculate gene means and std across all cells to be used for scaling
        gene_mean, gene_std = self.get_gene_mean_and_std(all_data)
        
        # Fit PCA embedding using randomized subsets of all_data
        self.log_trace('Fitting PCA model')
        pca_model = PCA(max_pca_components)
        for (i, cs) in enumerate(cell_subsets):
            # Scale features for a subset of cells and fit PCA
            all_data_subset = self.scale_features(csr_array(all_data[cs,]), gene_mean, gene_std)
            # Fit PCA model
            pca_model.fit(all_data_subset)
    
        # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
        x0, x1 = 0, min(len(pca_model.explained_variance_), max_pca_components - 1)
        y0, y1 = pca_model.explained_variance_[x0], pca_model.explained_variance_[x1]
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(pca_model.explained_variance_[:max_pca_components])]
        optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
        self.log_trace(f'  * Optimal PCA components: {optimal_pca_components}')
    
        # Elbow plot of explained variance
        plt.figure(figsize=(5, 3))
        plt.plot(range(x0 + 1, x1 + 2), pca_model.explained_variance_[:max_pca_components] / 100, 'b', marker='o', markersize=5, linewidth=1)
        plt.plot([optimal_pca_components, optimal_pca_components], [y0 / 100, y1 / 100], linestyle='dashed', color='red', linewidth=1)
        plt.xlabel('PC')
        plt.ylabel('Explained Variance (%)')
        plt.tight_layout()
        plt.savefig(f'{self.lower_dimensions_dir}/PCA_ExplainedVariance', dpi=300)
        plt.close()
        
        # Fit UMAP embedding using randomized subsets of all_data
        self.log_trace('Fitting UMAP model')
        umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
        for (i, cs) in enumerate(cell_subsets):
            # Scale features for a subset of cells
            all_data_subset = self.scale_features(csr_array(all_data[cs,]), gene_mean, gene_std)
            # PCA transform
            all_data_subset = pca_model.transform(all_data_subset)
            # Fit UMAP model
            umap_model.fit(all_data_subset[:, :optimal_pca_components])
        
        # Save gene means and std used for scaling (useful to add new cells if one does not want to refir PCA and UMAP)
        pd.DataFrame(gene_mean).to_csv(f'{self.lower_dimensions_dir}/gene_means.tsv.gz', sep='\t', header=False, index=False)
        pd.DataFrame(gene_std).to_csv(f'{self.lower_dimensions_dir}/gene_stds.tsv.gz', sep='\t', header=False, index=False)
        
        # Save PCA and UMAP models
        self.log_trace('Saving PCA and UMAP models to pickle file')
        pk.dump(pca_model, open(f'{self.lower_dimensions_dir}/pca_optimalcomp-{optimal_pca_components}.pkl', 'wb'))
        pk.dump(umap_model, open(f'{self.lower_dimensions_dir}/umap_neighbors-{neighbors}.pkl', 'wb'))
        
        # Save transformed data
        # N.B Scaled data is saved as a cell-to-gene full matrix, so it will be large
        self.log_trace('Saving embedded data')
        #scaled_out = gzip.open(f'{lower_dimensions_dir}/scaled_data.tsv.gz', 'wb')
        pca_out = gzip.open(f'{self.lower_dimensions_dir}/pca_data.tsv.gz', 'wb')
        umap_out = gzip.open(f'{self.lower_dimensions_dir}/umap_data.tsv.gz', 'wb')
        for i in range(all_data.shape[0] // self.max_cells_in_memory + 1):
            # Scale features for a subset of cells
            all_data_subset = self.scale_features(csr_array(all_data[i * self.max_cells_in_memory : (i + 1) * self.max_cells_in_memory,]), gene_mean, gene_std)
            out_text = '\n'.join(['\t'.join(line.astype(str)) for line in all_data_subset])
            out_text = out_text.encode() + (b'\n' if i < all_data.shape[0] // self.max_cells_in_memory else b'')
            #scaled_out.write(out_text)
            # PCA transform
            all_data_subset = pca_model.transform(all_data_subset)
            out_text = '\n'.join(['\t'.join(line.astype(str)) for line in all_data_subset])
            out_text = out_text.encode() + (b'\n' if i < all_data.shape[0] // self.max_cells_in_memory else b'')
            pca_out.write(out_text)
            # UMAP transform
            all_data_subset = umap_model.transform(all_data_subset[:, :optimal_pca_components])
            out_text = '\n'.join(['\t'.join(line.astype(str)) for line in all_data_subset])
            out_text = out_text.encode() + (b'\n' if i < all_data.shape[0] // self.max_cells_in_memory else b'')
            umap_out.write(out_text)
        
        # Close files
        #scaled_out.close()
        pca_out.close()
        umap_out.close()
        
        self.log_trace('All done!')
        self.log_trace(f'  * PCA and UMAP models were stored in {self.lower_dimensions_dir}')
        
    ### ------------------------------------ ###
    
    def scale_features(self, mtx: csr_array, gene_mean: np.ndarray=np.array([]), gene_std: np.ndarray=np.array([])) -> np.ndarray:
        
        # Z-score scaling of a sparse matrix
        # Note that after scaling the matrix will be dense since 0 values will be transformed to non-0
        
        # Get gene mean and std, if not provided
        if not len(gene_mean) or not len(gene_std):
            gene_mean, gene_std = self.get_gene_mean_and_std(mtx)
        
        # Scale features, parsing chunks of max_cells_in_memory
        mtx_scaled = (mtx.toarray() - gene_mean) / (gene_std + 1e-6)
        
        return mtx_scaled

    ### ------------------------------------ ###
    
    def umap_embedding(self, integrated_pca_path: str, neighbors: int=30, n_permutations: int=3) -> None:
        
        # N.B. PCA/UMAP models and embeddings are not kept in memory, but have to be loaded with self.load_reduced_dimensions afterwards
        self.log_trace('Reducing dataset dimensions')
        
        # Create dir to store dimensionality reduction data
        if not exists(self.lower_dimensions_dir):
            makedirs(self.lower_dimensions_dir)
        
        # Load integrated pca data
        try:
            pca_data = self.pca_data # From integrate_data function
        except:
            pca_data = pd.read_csv(integrated_pca_path, sep='\t', header=None).to_numpy()
        
        # Create subsets of self.max_cells_in_memory cells to fit the PCA and UMAP models to
        # N.B. All cells are passed through 3 times in different groupings
        if pca_data.shape[0] < self.max_cells_in_memory:
            cell_subsets = self.make_cell_subsets(pca_data.shape[0], self.max_cells_in_memory)
        else:
            cell_subsets = [subset for _ in range(n_permutations) for subset in self.make_cell_subsets(pca_data.shape[0], self.max_cells_in_memory)]
        
        # Fit UMAP embedding using randomized subsets of all_data
        self.log_trace('Fitting UMAP model')
        umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
        for (i, cs) in enumerate(cell_subsets):
            # Subset cells
            all_data_subset = pca_data[cs,]
            # Fit UMAP model
            umap_model.fit(all_data_subset)
        
        # Save UMAP model
        self.log_trace('Saving UMAP model to pickle file')
        pk.dump(umap_model, open(f'{self.lower_dimensions_dir}/umap_neighbors-{neighbors}.pkl', 'wb'))
        
        # Save transformed data
        pd.DataFrame(pca_data).to_csv(f'{self.lower_dimensions_dir}/pca_data.tsv.gz', sep='\t', header=False, index=False)
        
        self.log_trace('Saving embedded data')
        umap_out = gzip.open(f'{self.lower_dimensions_dir}/umap_data.tsv.gz', 'wb')
        for i in range(pca_data.shape[0] // self.max_cells_in_memory + 1):
            # UMAP transform
            umap_data_subset = umap_model.transform(pca_data[i * self.max_cells_in_memory : (i + 1) * self.max_cells_in_memory,])
            out_text = '\n'.join(['\t'.join(line.astype(str)) for line in umap_data_subset])
            out_text = out_text.encode() + (b'\n' if i < self.all_data.shape[0] // self.max_cells_in_memory else b'')
            umap_out.write(out_text)
        
        # Close files
        umap_out.close()

        self.log_trace('All done!')
        self.log_trace(f'  * PCA and UMAP models were stored in {self.lower_dimensions_dir}', False)
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_gene_mean_and_std(mtx: csr_array) -> tuple[np.ndarray, np.ndarray]:
        
        # Transposing matrix from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_array slicing is efficient along rows
        mtx_t = mtx.transpose()
        
        # Calculate features mean and variance
        # Remember that variance = mean of squared values - square of mean values
        mtx_squared = mtx_t.copy()
        mtx_squared.data **= 2
        gene_mean = np.asarray(mtx_t.mean(axis=1)).reshape(-1)
        gene_std = (np.asarray(mtx_squared.mean(axis=1)).reshape(-1) - (gene_mean ** 2))**0.5
        
        # Reshaping mean and std
        gene_mean = gene_mean.reshape((1, gene_mean.shape[0]))
        gene_std = gene_std.reshape((1, gene_std.shape[0]))
        
        return gene_mean, gene_std
    
    ### ------------------------------------ ###

    @staticmethod
    def make_cell_subsets(N: int, max_cells: int) -> list[np.ndarray]:
        
        subsets, indexes = [], np.arange(N)
        while len(indexes):
            new_subset = np.random.choice(indexes, min(max_cells, len(indexes)), replace=False)
            subsets.append(new_subset)
            indexes = indexes[~ np.isin(indexes, new_subset)]
        
        return subsets
    
    ### ------------------------------------ ###
    ### CLUSTERING                           ###
    ### ------------------------------------ ###

    def check_cluster_composition(self) -> None:
        
        try:
            _ = self.clusters
        except:
            self.log_trace('ERROR: need to run/load clustering analysis first')
        
        # Get unique clusters
        clusters_id = np.sort(np.unique(self.clusters))
        
        # Get unique datasets
        datasets_labels = self.datasets_metadata.datasets_labels.to_list()
        
        # Get total number of cells per dataset
        datasets_cells = self.datasets_metadata.cell_num.to_list()
        
        # Init dict for results
        results = {key : [] for key in ['cluster', 'cluster_cells'] + datasets_labels + ['hypergeom_test_flag']}
        results['cluster'] = list(clusters_id)
        
        # Check number of cells of each dataset in individual clusters
        for cl in clusters_id:
            # Get total number of cells in cluster
            cl_cells = sum(self.clusters == cl)
            results['cluster_cells'].append(cl_cells)
            # Check cluster composition
            pval_flag = 0
            for ds_l, ds_c in zip(datasets_labels, datasets_cells):
                # Number of cells from dataset ds_l in cluster cl
                ds_cells = sum([1 for cell in np.array(self.all_cells)[self.clusters == cl] if cell.endswith(f'_{ds_l}')])
                # Hypergeometric test
                pval = hypergeom(sum(datasets_cells), cl_cells, ds_c).sf(ds_cells - 1)
                if pval < 0.05:
                    pval_flag = 1
                # Store results
                results[ds_l].append(ds_cells / cl_cells)
            results['hypergeom_test_flag'].append(pval_flag)
            
        # Convert results to pandas dataframe
        results = pd.DataFrame(results)
        
        # Save results
        results.to_csv(f'{self.clusters_dir}/cluster_composition.tsv', sep='\t', header=True, index=False)
    
    ### ------------------------------------ ###
        
    def cluster_cells(self, n_neighbors: int=10, resolution_parameter: float=1.0, beta: float=0.01, run_optimizer: bool=True) -> None:
        
        self.log_trace('Clustering cells')
        
        # Create dir to store cluster embeddings
        if not exists(self.clusters_dir):
            makedirs(self.clusters_dir)
        
        # Setting random seed (helps with clustering consistency)
        random.seed(42)
        
        # Check that self.lower_dimensions_dir exists
        try:
            _ = self.lower_dimensions_dir
        except:
            self.log_trace('  * ERROR: missing path to lower dimension data (self.lower_dimensions_dir)', False)
            return None
        
        # Load PCA and UMAP embdeddings
        _, pca_data, _, _, _ = self.load_reduced_dimensions(self.lower_dimensions_dir)

        # Run optimizer if run_optimizer is True
        if run_optimizer:
            optimizer_instance = LeidenClusteringEvaluator()
            optimizer_params = {
                'n_neighbors_vals' : [10, 20, 30, 40, 50],
                'resolution_vals' : [0.4, 0.6, 0.8, 1.0, 1.2],
                'subset_fraction' : 0.75,
                'n_iterations' : 10
            }
            evaluation = optimizer_instance.evaluate_parameters(pd.DataFrame(pca_data, columns=[f'PCA{i+1}' for i in range(pca_data.shape[1])]), **optimizer_params)
            evaluation.to_csv(f'{self.clusters_dir}/clustering_parameter_optimization.tsv', sep='\t', header=True, index=False)
            optimal_clustering_params = evaluation.sort_values('ari', ascending=False).iloc[0]
            n_neighbors, resolution_parameter = optimal_clustering_params[['n_neighbors', 'resolution']]
            n_neighbors = int(n_neighbors)
            self.log_trace('Optimal parameters found:')
            self.log_trace(f'n_neighbors: {n_neighbors}')
            self.log_trace(f'resolution: {resolution_parameter}')

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
        clusters = graph.community_leiden(objective_function='modularity', weights=None, resolution_parameter=resolution_parameter, beta=beta, initial_membership=None, n_iterations=2, node_weights=None).membership
        clusters = np.array(clusters)
        
        self.clusters = clusters
        
        # Save clusters to file
        pd.DataFrame(clusters).to_csv(f'{self.clusters_dir}/{self.cluster_embeddings_file}', sep='\t', header=False, index=False)
        
        # Check cluster composition
        self.check_cluster_composition()
        
        self.log_trace('All done!')
        self.log_trace(f'  * Clusters info was stored in {self.clusters_dir}/{self.cluster_embeddings_file}', False)
     
    ### ------------------------------------ ###
    
    def find_clusters_markers(self, min_fc: float=0.25, min_pct: float=0.1, min_cells_in_cluster: int=50) -> None:

        self.log_trace('Finding cluster markers')

        # Create dir to store cluster markers
        if not exists(self.clusters_dir):
            makedirs(self.clusters_dir)
        
        # Init results table
        cluster_markers = {"cluster" : [], "gene_name" : [], "cluster_pct" : [], "other_clusters_pct" : [], "cluster_expression" : [], "other_clusters_expression" : [], "log2FC" : [], "pval" : [], "padj" : []}
        
        # Get clusters ids
        clusters_id = np.sort(np.unique(self.clusters))
        
        # Find markers in each cluster vs other clusters
        for cl in clusters_id:
            # Extracting cells indices
            cluster_cells = np.where(self.clusters == cl)[0]
            other_cells = np.where(self.clusters != cl)[0]
            if len(cluster_cells) < min_cells_in_cluster: # Skip clusters with too few cells
                self.log_trace(f'  * WARNING: cluster {cl} comprises only {len(cluster_cells)} cells. Skipping cluster...', False)
                continue
            self.log_trace(f'  * Finding markers of cluster {cl}')
            # Running analysis
            significant_genes_name, log2_fc, pop1_expression, pop2_expression, pct_pop1, pct_pop2, pvals, padjs = self.compare_populations(self.all_genes, self.all_data, cluster_cells, other_cells, self.max_cells_in_memory, min_fc, min_pct)
            # Storing data
            cluster_markers["cluster"].extend([cl for _ in range(len(significant_genes_name))])
            cluster_markers["gene_name"].extend(significant_genes_name)
            cluster_markers["cluster_pct"].extend(pct_pop1)
            cluster_markers["other_clusters_pct"].extend(pct_pop2)
            cluster_markers["cluster_expression"].extend(pop1_expression)
            cluster_markers["other_clusters_expression"].extend(pop2_expression)
            cluster_markers["log2FC"].extend(log2_fc)
            cluster_markers["pval"].extend(pvals)
            cluster_markers["padj"].extend(padjs)
        cluster_markers = pd.DataFrame(cluster_markers)
        self.cluster_markers = cluster_markers
        
        # Save to file
        cluster_markers.to_csv(f'{self.clusters_dir}/{self.cluster_markers_file}', sep='\t', header=True, index=False)
        
        self.log_trace('All done!')
        self.log_trace(f'  * Clusters markers were saved to {self.clusters_dir}/{self.cluster_markers_file}', False)
    
    ### ------------------------------------ ###
    
    def score_cluster_identity_cipr(self, ref_expression: pd.DataFrame, ref_samples: pd.DataFrame, clusters_id: list[int]=[]) -> None:
        
        # Similar to CIPR, but simpler and better
        self.log_trace('Assigning identities to clusters')
        
        # Create dir to store dimensionality reduction data
        if not exists(self.clusters_dir):
            makedirs(self.clusters_dir)
        
        # Get clusters ids
        if not len(clusters_id):    
            clusters_id = np.sort(np.unique(self.clusters)).tolist()

        # Listing possible cell types and subtypes
        cell_classes = ref_samples['reference_cell_type'].toarray()
        cell_subclasses = ref_samples['long_name'].toarray()
        
        # Init cluster identity matrix
        cluster_identities = {'cluster' : [],
                              'cell_type' : [],
                              'cell_type_prob' : [],
                              'cell_sub_type' : [],
                              'cell_sub_type_prob' : []}
        
        # Assessing identity probabilities for each cluster
        for cl in clusters_id:
            self.log_trace(f'Scoring cluster {cl} identity')
            # Subset markers
            cl_markers = self.cluster_markers.loc[
                (self.cluster_markers.cluster == cl) &
                (self.cluster_markers.padj < 0.05) &
                (self.cluster_markers.gene_name.str.lower().isin(ref_expression.iloc[:, 0].str.lower()))
            ]
            cl_markers = cl_markers.loc[cl_markers.gene_name.str.lower().duplicated() == False,]
            cl_markers = cl_markers.iloc[np.argsort(cl_markers.gene_name.str.lower()),]
            
            # Testing cell_types
            probs = []
            for col in ref_expression.columns[1:]:
                ref_logfc = ref_expression.loc[ref_expression.gene_name.str.lower().isin(cl_markers.gene_name.str.lower()), ["gene_name", col]]
                ref_logfc = ref_logfc.iloc[np.argsort(ref_logfc.gene_name.str.lower()),]
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
                class_prob = probs[ref_samples['reference_cell_type'] == c_class].sum()
                if class_prob > best_class_prob:
                    best_class = c_class
                    best_class_prob = class_prob

            # Finding most probable subclass
            subclass_probs = probs[ref_samples['reference_cell_type'] == best_class]
            best_subclass_prob = subclass_probs.max()
            best_subclass = list(cell_subclasses[ref_samples['reference_cell_type'].values == best_class][subclass_probs == best_subclass_prob])[0]
            
            # Adding cluster info
            cluster_identities['cluster'].append(cl)
            cluster_identities['cell_type'].append(best_class)
            cluster_identities['cell_type_prob'].append(best_class_prob)
            cluster_identities['cell_sub_type'].append(best_subclass)
            cluster_identities['cell_sub_type_prob'].append(best_subclass_prob)
        cluster_identities = pd.DataFrame(cluster_identities)
        self.cluster_identities = cluster_identities
        
        # Saving to file
        cluster_identities.to_csv(f'{self.clusters_dir}/{self.cluster_identities_file}', sep='\t', header=True, index=False)
        
        self.log_trace('All done!')
        self.log_trace(f'  * Clusters identities were saved to {self.clusters_dir}/{self.cluster_identities_file}', False)

    ### ------------------------------------ ###

    def score_cluster_identity_sctype(self, markers: pd.DataFrame, tissues: list[str], clusters_id: list[int]=[]) -> None:

        # Python implementation of ScType
        self.log_trace('Assigning identities to clusters')
        
        # Create dir to store dimensionality reduction data
        if not exists(self.clusters_dir):
            makedirs(self.clusters_dir)
        
        # Get clusters ids
        if not len(clusters_id):    
            clusters_id = np.sort(np.unique(self.clusters)).tolist()

        # Init sc-type
        sc_type_instance = ScTypeLike(self.all_data, self.all_genes, markers, tissues)
        self.log_trace('  * Scoring individual cells')
        sc_type_instance.score_cell_identity()
        self.log_trace('  * Collaping to cluster level')
        sc_type_instance.reduce_to_cluster_level(self.clusters.tolist())
        cluster_identities = sc_type_instance.cluster_identities
        cluster_identities.columns = ['cluster', 'cell_type']
        self.cluster_identities = cluster_identities
        
        # Saving to file
        cluster_identities.to_csv(f'{self.clusters_dir}/{self.cluster_identities_file}', sep='\t', header=True, index=False)
        
        self.log_trace('All done!')
        self.log_trace(f'  * Clusters identities were saved to {self.clusters_dir}/{self.cluster_identities_file}', False)

    ### ------------------------------------ ###
    
    @staticmethod
    def compare_populations(gene_names: list[str], data: csr_array, pop1: np.ndarray, pop2: np.ndarray, max_cells_in_memory: int=10000, min_fc: float=0.25, min_pct: float=0.1) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        
        # Split all_data into pop1 and pop2
        pop1_data = csr_array(data[pop1, ])
        pop2_data = csr_array(data[pop2, ])
        
        # Transposing matrices from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_array slicing is efficient along rows
        pop1_data = pop1_data.transpose()
        pop2_data = pop2_data.transpose()
        
        # Extracting average exression values and log2Fc
        pop1_expression = np.asarray(pop1_data.mean(axis=1)).reshape(-1)
        pop2_expression = np.asarray(pop2_data.mean(axis=1)).reshape(-1)
        log2_fc = np.log2(pop1_expression / (pop2_expression + 0.001) + 0.001)
        
        # Calculating the number of cells in which a certain feature is detected (i.e. normalized_count > 0)
        pct_pop1 = np.diff(pop1_data.T.tocsc().indptr) / pop1_data.shape[1]
        pct_pop2 = np.diff(pop2_data.T.tocsc().indptr) / pop2_data.shape[1]
        
        # Filtering features to test by min_fc and min_pct
        features_to_test = (abs(log2_fc) > min_fc) & ((pct_pop1 > min_pct) | (pct_pop2 > min_pct))
        features_to_test = np.where(features_to_test)[0]
        
        # Subset and transpose matrices
        pop1_data = pop1_data.T[:, features_to_test]
        pop2_data = pop2_data.T[:, features_to_test]
        
        # Significance testing and BH correction
        if pop1_data.shape[0] < max_cells_in_memory and pop2_data.shape[0] < max_cells_in_memory:
            pvals = mannwhitneyu(pop1_data.toarray(), pop2_data.toarray(), alternative="two-sided", axis=0).pvalue
        else:
            pvals = []
            for gene in range(pop1_data.shape[1]):
                pval = mannwhitneyu(pop1_data[:, gene].toarray(), pop2_data[:, gene].toarray(), alternative="two-sided", axis=0).pvalue[0]
                pvals.append(pval)
            pvals = np.array(pvals)
        _, padjs = fdrcorrection(pvals, alpha=0.05, is_sorted=False)
        
        # Filtering
        significant_genes_name = np.array(gene_names)[features_to_test][padjs < 0.05].tolist()
        log2_fc = log2_fc[features_to_test][padjs < 0.05]
        pop1_expression = pop1_expression[features_to_test][padjs < 0.05]
        pop2_expression = pop2_expression[features_to_test][padjs < 0.05]
        pct_pop1 = pct_pop1[features_to_test][padjs < 0.05]
        pct_pop2 = pct_pop2[features_to_test][padjs < 0.05]
        pvals = pvals[padjs < 0.05]
        padjs = padjs[padjs < 0.05]
        
        return significant_genes_name, log2_fc, pop1_expression, pop2_expression, pct_pop1, pct_pop2, pvals, padjs

    ### ------------------------------------ ###
    ### TRAJECTORY ANALYSIS                  ###
    ### ------------------------------------ ###
    
    def create_mst(self, clusters_id: list[int]=[]) -> None:
        
        # Kruskal Minimal Spanning Tree algorithm
        self.log_trace('Starting trajectory analysis')
        try:
            _ = self.lower_dimensions_dir
        except:
            self.log_trace('  * ERROR: missing path to lower dimension data (self.lower_dimensions_dir)', False)
            return None
        
        # Create dir to store dimensionality reduction data
        if not exists(self.trajectories_dir):
            makedirs(self.trajectories_dir)
        
        # Load PCA embdeddings
        _, pca_data, _, _, _ = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Finding clusters centroids (graph vertices)
        if not len(clusters_id):
            clusters_id = np.sort(np.unique(self.clusters)).tolist()
        centroids = np.array([pca_data[self.clusters == cl, ].mean(axis = 0) for cl in clusters_id])
        
        # Calculating pairwise distances between cluster centroids (graph edges)
        distances = []
        for c1_index,c1 in enumerate(centroids):
            for c2_index,c2 in enumerate(centroids[c1_index+1:]):
                new_distance = (((c1 - c2)**2).sum())**0.5
                distances.append([c1_index, c1_index + c2_index + 1, new_distance])
        distances = np.array(distances)
        distances = distances[distances[:, 2].argsort(),]
        
        # Init parents and ranks
        n_vertices = len(clusters_id)
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
            branches[br_n] = [clusters_id[index] for index in br]
        self.branches = branches
        
        # Save branches to file
        out_text = '\n'.join(['\t'.join([str(b) for b in br]) for br in branches])
        with open(f'{self.trajectories_dir}/{self.trajectories_file}', 'w') as output:
            output.write(out_text)
        
        self.log_trace('All done!')
        self.log_trace(f'  * Trajectory analysis was saved to {self.trajectories_dir}', False)

    ### ------------------------------------ ###
    
    def generate_pseudotime(self) -> None:
        
        self.log_trace('Calculating cells pseudotime along trajectories')
        try:
            _ = self.lower_dimensions_dir
        except:
            self.log_trace('  * ERROR: missing path to lower dimension data (self.lower_dimensions_dir)', False)
            return None
        try:
            _ = self.trajectories_dir
            _ = self.trajectory_file_name
        except:
            self.log_trace('. * ERROR: please run or load trajectory analysis', False)
            return None
        
        # Load PCA embdeddings
        pca_model, pca_data, _, _, _ = self.load_reduced_dimensions(self.lower_dimensions_dir)
        pca_explained_variance = pca_model.explained_variance_
        pca_explained_variance = (pca_explained_variance * 100 / pca_explained_variance.sum())[:pca_data.shape[1]]
        
        # Listing clusters used
        clusters_id = list({vertex for br in self.branches for vertex in br})
        clusters_id.sort()
        
        # Calculating cluster centroids
        centroids = np.array([pca_data[self.clusters == cl, ].mean(axis=0) for cl in clusters_id])
        
        # Calculating pseudotime for each cell
        pseudotimes = []
        for c in range(pca_data.shape[0]):
            coords = pca_data[c, ]
            cluster = self.clusters[c]
            branch = [br for br in self.branches if cluster in br]
            if not len(branch) or cluster not in clusters_id: # Cluster is not connected
                pseudotimes.append(0)
                continue
            branch = branch[0]
            vertex_index = branch.index(cluster)
            if vertex_index == 0: # Cell belongs to root vertex. Pseudotime = projection on first edge of branch
                # Projecting cell to first edge of branch
                branch_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[:2]])
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
                # Computing pseudotime
                #cell_pseudotime = ((projection**2).sum())**0.5
                cell_pseudotime = (((projection * pca_explained_variance)**2).sum())**0.5
            elif vertex_index == len(branch) - 1: # Cell belongs to end vertex. Pseudotime = projection on last edge of branch
                # Projecting cell to last edge of branch
                branch_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[-2:]])
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
                # Calculating previous branch length
                previous_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[:vertex_index]])
                #added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
                added_distance = sum((((np.diff(previous_centroids) * pca_explained_variance)**2).sum(axis=1))**0.5)
                # Computing pseudotime
                cell_pseudotime = (((projection * pca_explained_variance)**2).sum())**0.5 + added_distance
            else: # Cell belongs to a end/middle vertex. Pseudotime = projection on edge of branch + previous branches' length
                # Finding closest cluster centroid to the cell that the cell doesn't belong to
                near_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[vertex_index - 1 : vertex_index + 3]])
                distance_1 = (((near_centroids[0] - coords)**2).sum())**0.5
                distance_2 = (((near_centroids[-1] - coords)**2).sum())**0.5
                # Projecting cell to nearest edge of branch
                branch_centroids = near_centroids[:2] if distance_1 < distance_2 else near_centroids[-2:]
                vertex_index = vertex_index if distance_1 < distance_2 else vertex_index + 1
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
                # Calculating previous branch length
                previous_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[:vertex_index]])
                #added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
                added_distance = sum((((np.diff(previous_centroids) * pca_explained_variance)**2).sum(axis=1))**0.5)
                # Computing pseudotime
                cell_pseudotime = (((projection * pca_explained_variance)**2).sum())**0.5 + added_distance
            pseudotimes.append(cell_pseudotime)
        
        # Normalize pseudotime to 0-100
        pseudotimes = (np.array(pseudotimes) - min(pseudotimes)) * 100 / (max(pseudotimes) - min(pseudotimes))
        self.pseudotimes = pseudotimes
      
        # Save pseudotime to path
        out_path = f'{self.trajectories_dir}/{self.pseudotime_file}'
        pd.DataFrame(pseudotimes).to_csv(out_path, sep='\t', header=False, index=False)
        
        self.log_trace('All done!')
        self.log_trace(f'  * Trajectory analysis was saved to {self.trajectories_dir}', False)
    
    ### ------------------------------------ ###
    ### CELL CYCLE ANALYSIS                  ###
    ### ------------------------------------ ###
    
    def score_cell_cycle(self, g1pm_features: list[str], s_features: list[str], g2m_features: list[str], gene_pool: list[str]=[], bins: int=50, ctrl_genes_num: int=50) -> None:
        
        self.log_trace('Scoring cell cycle')
        
        # Create dir to store dimensionality reduction data
        if not exists(self.cell_cycle_dir):
            makedirs(self.cell_cycle_dir)

        # Define cell cycle phases
        phases = ["G1PM", "S", "G2M"]
        
        # Scoring gene sets
        toggle = 1
        cell_cycle_scores = []
        for set_name, gene_set in zip(phases, [g1pm_features, s_features, g2m_features]):
            self.log_trace(f'  * Scoring {set_name} phase', False)
            set_score = self.score_gene_set(gene_set, gene_pool, bins, ctrl_genes_num)
            if set_score is None:
                toggle = 0
                bad_gene_set = set_name
                break
            else:
                cell_cycle_scores.append(set_score)
        if toggle == 0:
            self.log_trace(f'Error: {bad_gene_set} list is bad')
            return None
        else:
            cell_cycle_scores = np.stack(cell_cycle_scores).T
            cell_cycle_scores = pd.DataFrame(cell_cycle_scores, columns=[f'{p}_score' for p in phases])
            cell_cycle_scores["cell_cycle_phase"] = [phases[x] if cell_cycle_scores.values[y, x] > 0 else "G1" for y,x in enumerate(np.argmax(cell_cycle_scores.to_numpy(), axis=1))]
            self.cell_cycle_scores = cell_cycle_scores
            # Save cell cycle scores to file
            out_path = f'{self.cell_cycle_dir}/{self.cell_cycle_scores_file}'
            cell_cycle_scores.to_csv(out_path, sep='\t', header=True, index=False)

            self.log_trace('All done!')
            self.log_trace(f'  * Cell cycle scores were saved to {out_path}', False)
    
    ### ------------------------------------ ###

    def score_gene_set(self, gene_set: list[str], gene_pool: list[str]=[], bins: int=50, ctrl_genes_num: int=50):

        # Subsetting gene_set and gene_list for detected genes
        gene_set = [gene for gene in gene_set if gene in self.all_genes]
        gene_pool = [gene for gene in gene_pool if gene in self.all_genes]

        # Making sure that gene_set and gene_pool have len > 0
        if not len(gene_set):
            self.log_trace("Error: empty gene_set")
            return None

        # Checking gene pool
        if not len(gene_pool):
            gene_pool = self.all_genes

        # Compute average gene expression levels
        gene_means = np.array(self.all_data.mean(axis=0)).ravel()
        
        # Rank genes based on binned expression level, then for each bin of the genes in the gene_set, pick ctrl_genes random genes for genes with matched binned expression
        bin_size = len(gene_means) // bins
        expression_order = np.argsort(np.argsort(gene_means))
        set_ranks = expression_order[np.isin(self.all_genes, gene_set)].tolist()
        ctrl_indexes = []
        for rank in set_ranks:
            random_pool = np.where((expression_order >= (rank - bin_size // 2)) &
                                   (expression_order <= (rank + bin_size // 2)) &
                                   ~(np.isin(expression_order, set_ranks)))[0].tolist()
            np.random.shuffle(random_pool)
            ctrl_indexes.extend(random_pool[:ctrl_genes_num])
        ctrl_indexes = list(set(ctrl_indexes)) # Removing duplicates

        # Computing the mean of scaled expression values of gene_set genes for each cell
        set_indexes = [self.all_genes.index(g) for g in gene_set if g in self.all_genes]
        set_means = self.scale_features(csr_array(self.all_data[:, set_indexes])).mean(axis=1)
        #if self.all_data.shape[0] <= self.max_cells_in_memory:
        #    set_means = self.scale_features(csr_array(self.all_data[:, set_indexes])).mean(axis=1)
        #else:
        #    set_means = np.array([self.scale_features(csr_array(self.all_data[:, g])).reshape(-1) for g in set_indexes]).T.mean(axis=1)

        # Computing the mean of scaled expression values of  ctrl_set genes for each cell
        ctrl_means = self.scale_features(csr_array(self.all_data[:, ctrl_indexes])).mean(axis=1)
        #if self.all_data.shape[0] <= self.max_cells_in_memory:
        #    ctrl_means = self.scale_features(csr_array(self.all_data[:, ctrl_indexes])).mean(axis=1)
        #else:
        #    ctrl_means = np.array([self.scale_features(csr_array(self.all_data[:, g])).reshape(-1) for g in ctrl_indexes]).T.mean(axis=1)
        
        set_score = set_means - ctrl_means
        
        return set_score

    ### ------------------------------------ ###
    ### PSEUDOBULK                           ###
    ### ------------------------------------ ###

    def make_pseudobulk(self, cells: list[int], method: str='median') -> pd.Series:

        data_sub = csr_array(self.all_data[cells]).toarray()
        if method == 'median':
            pseudobulk = pd.Series(np.median(data_sub, axis=0), index=self.all_genes)
        elif method == 'mean':
            pseudobulk = pd.Series(np.mean(data_sub, axis=0), index=self.all_genes)
        else:
            self.log_trace(f'ERROR: Unknown method for pseudobulk generation "{method}"')
            return pd.Series()

        return pseudobulk

    ### ------------------------------------ ###
    ### LIGAND RECEPTOR INTERACTIONS         ###
    ### ------------------------------------ ###

    def analyze_cell_interactions(self, interaction_pairs: pd.DataFrame, cells: list[int]=[], iterations: int=10000) -> None:

        self.log_trace('Running cell interactions analysis')
        self.log_trace('Prepping data')

        # Create dir to store sparse matrices after preprocessing
        if not exists(self.cell_interactions_dir):
            makedirs(self.cell_interactions_dir)

        # Label each cell with identities
        cl_ti_ct = {cl : ct for _,(cl,ct) in self.cluster_identities.iterrows()}
        cell_identities = np.array([cl_ti_ct[cl] for cl in self.clusters])

        # Subset count data rows (cells) and update list of identities
        if len(cells):
            data_sub = csr_array(self.all_data[cells])
            cell_identities = cell_identities[cells]
        else:
            data_sub = self.all_data.copy()

        # Filter for detected genes
        interaction_pairs = interaction_pairs.loc[
            (interaction_pairs['receptor'].isin(self.all_genes)) &
            (interaction_pairs['ligand'].isin(self.all_genes))
        ]
        
        # Subset count data columns (genes) to make calculations faster
        unique_ligands_and_receptors = np.unique(interaction_pairs.to_numpy().ravel()).tolist()
        data_sub = csr_array(data_sub[:, [self.all_genes.index(lr) for lr in unique_ligands_and_receptors]])

        # Score observed interactions scores
        self.log_trace('Calculating observed interaction scores')
        obs_scores = self.score_cell_interactions(data_sub, cell_identities, interaction_pairs, unique_ligands_and_receptors)

        # Initialize dictionary to accumulate permutation scores for the null background
        shuffle_score_lists = {key: [] for key in obs_scores.keys()}

        # Run iterations to approximate null distribution
        self.log_trace(f'Running {iterations} shuffling permutations for Z-test')
        cell_identities_shuffle = cell_identities.copy()
        for I in range(iterations):
            if (I + 1) % 100 == 0 and I > 0:
                print(f'Ran {I+1} iterations', end='\r')
            # Shuffle identity labels while preserving cluster sizes
            np.random.shuffle(cell_identities_shuffle)
            # Score interactions after shuffling
            shuffle_scores = self.score_cell_interactions(data_sub, cell_identities_shuffle, interaction_pairs, unique_ligands_and_receptors)
            # Append to null distribution lists
            for key in obs_scores.keys():
                shuffle_score_lists[key].append(shuffle_scores.get(key, 0.0))

        # Calculate Z-scores and significance
        self.log_trace('Running statistics')
        cell_interactions, cell_interactions_header = [], ['sender', 'receiver', 'ligand', 'receptor', 'observed_score', 'mean_shuffle_score', 'z_score', 'pval']
        for key,s_obs in obs_scores.items():
            # Extract details of ligand/receptor pair and find related null distribution
            receptor, ligand, receiver, sender = key
            null_distribution = np.array(shuffle_score_lists[key])
            # Null distribution stats
            mean_null = np.mean(null_distribution)
            std_null = np.std(null_distribution)
            # Z-score (safe)
            if std_null > 0:
                z_score = (s_obs - mean_null) / std_null
            else:
                z_score = 0.0
            # P-value estimation (empirical) (faster than scipy.stats)
            pval = np.sum(null_distribution >= s_obs) / iterations
            cell_interactions.append((sender, receiver, ligand, receptor, s_obs, mean_null, z_score, pval))
        cell_interactions = pd.DataFrame(cell_interactions, columns=cell_interactions_header)

        # Benjamini-Hochberg adjustment
        self.log_trace('Adjusting p-values')
        cell_interactions.loc[:, 'padj'] = fdrcorrection(cell_interactions['pval'], alpha=0.05, is_sorted=False)[1]
        self.cell_interactions = cell_interactions
        cell_interactions.to_csv(f'{self.cell_interactions_dir}/{self.cell_interactions_file}', sep='\t', index=False)

        self.log_trace('All done!')
        self.log_trace(f'  * Cell interactions scores were saved to {self.cell_interactions_dir}/{self.cell_interactions_file}', False)

    ### ------------------------------------ ###

    @staticmethod
    def score_cell_interactions(data: csr_array, cell_identities: NDArray, receptor_ligand_pairs: pd.DataFrame, unique_ligands_and_receptors: list[str]=[]) -> dict[tuple[str, str, str, str], float]:

        # Extract unique cell types and ligands/receptors
        if not len(unique_ligands_and_receptors):
            unique_ligands_and_receptors = np.unique(receptor_ligand_pairs.to_numpy().ravel()).tolist()
        unique_cell_identities = list(set(cell_identities))

        # Init pd.DataFrame of cell types' mean normalized expressions
        cell_type_norm_means = pd.DataFrame(
            np.zeros((len(unique_cell_identities), len(unique_ligands_and_receptors)), dtype=float),
            index=unique_cell_identities,
            columns=unique_ligands_and_receptors
        )

        # Compute mean expression normalized with detection fraction for all cell types
        for ct in unique_cell_identities:
            # Find cells indices and subset counts
            ct_idxs = np.where(cell_identities == ct)[0]
            ct_data = csr_array(data[ct_idxs,])
            # Compute mean normalized by detection fraction
            ct_norm_means = ct_data.mean(axis=0) * (ct_data > 0).mean(axis=0)
            cell_type_norm_means.loc[ct, :] = ct_norm_means
        
        # Compute scores
        scores = {}
        for _,(r,l) in receptor_ligand_pairs.iterrows():
            for ct1 in unique_cell_identities:
                for ct2 in unique_cell_identities:
                    score = cell_type_norm_means.loc[ct1, r] * cell_type_norm_means.loc[ct2, l]
                    scores[(r, l, ct1, ct2)] = score

        return scores

    ### ------------------------------------ ###
    ### PLOTS                                ###
    ### ------------------------------------ ###
    
    def check_plot_dir(self) -> None:

        # Sanity check
        if not hasattr(self, 'plot_dir'):
            toggle = 1
        elif not exists(self.plot_dir) or not isdir(self.plot_dir):
            toggle = 1
        else:
            toggle = 0

        if toggle == 1:    
            # Create dir to store plots
            self.plot_dir = f'{self.repo_name}/data_plots'
            if not exists(self.plot_dir):
                makedirs(self.plot_dir)
    
    ### ------------------------------------ ###
    
    def plot_cell_cycle(self, datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters) or not hasattr(self, 'clusters'):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Adding expression data, then sorting by smallest value
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['CellCyclePhase'] = self.cell_cycle_scores.loc[:, "cell_cycle_phase"].to_list()
        
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
        
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='CellCyclePhase', hue_order=['G1PM', 'G1', 'S', 'G2M'], palette='Set1', marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Cell Cycle Phase', markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{self.plot_dir}/Cell_Cycle_Phase_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_cell_type(self, use_subtype: bool=False, cell_types: list[str]=[], datasets: list[str]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not use_subtype:
            cl_to_id = {cl : cl_id for _,(cl,cl_id) in self.cluster_identities[['cluster', 'cell_type']].iterrows()}
            cell_identities = np.array([cl_to_id[cl] for cl in self.clusters])
            plot_name = f'{self.plot_dir}/Cell_type_UMAP.png'
        else:
            cl_to_id = {cl : cl_id for _,(cl,cl_id) in self.cluster_identities[['cluster', 'cell_sub_type']].iterrows()}
            cell_identities = np.array([cl_to_id[cl] for cl in self.clusters])
            plot_name = f'{self.plot_dir}/Cell_subtype_UMAP.png'
        if not len(cell_types):
            identity_filter = np.array([True for _ in range(self.all_data.shape[0])])
            cell_types = np.sort(np.unique(cell_identities)).tolist()
        else:
            identity_filter = np.isin(cell_identities, cell_types)
        cell_filter = list(dataset_filter & identity_filter)
        
        # Subsetting cluster color palette
        color_palette = self.cluster_colors[:len(cell_types)]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(cell_types) / 16) + 1
        
        # Adding cluster data, then sorting by smallest value
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['Cell type'] = cell_identities

        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]

        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Cell type', hue_order=cell_types, palette=color_palette, marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Cell type', ncol=legend_cols, markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(plot_name, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_clusters(self, datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
            clusters = np.sort(np.unique(self.clusters)).tolist()
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Subsetting cluster color palette
        color_palette = [self.cluster_colors[cl] for cl in clusters]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(clusters) / 16) + 1
        
        # Adding cluster data, then sorting by smallest value
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['Clusters'] = self.clusters.tolist()
    
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
    
        # Finding cluster centers for writing text
        centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
    
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Clusters', palette=color_palette, marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Clusters', ncol=legend_cols, markerscale=markerscale)
        for c_num,cl in enumerate(centers):
            x, y = cl
            plt.text(x, y, str(clusters[c_num]), horizontalalignment='center', size='small', color='black', weight='semibold')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{self.plot_dir}/Clusters_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_datasets(self, datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
            datasets = self.datasets_metadata.datasets_labels.to_list()
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters) or not hasattr(self, 'clusters'):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Subsetting cluster color palette
        color_palette = self.cluster_colors[:len(datasets)]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(datasets) / 16) + 1
        
        # Adding dataset data, then sorting by smallest value
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['Datasets'] = [cell.split('_')[-1] for cell in self.all_cells]
    
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
    
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Datasets', palette=color_palette, marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Datasets', ncol=legend_cols, markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{self.plot_dir}/Datasets_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_gene_expression(self, target_gene: str, datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Extracting expression data for target_gene
        if target_gene not in list(self.all_genes):
            toggle = False
        else:
            toggle = True
            gene_index = self.all_genes.index(target_gene)
            expression = csr_array(self.all_data[:, [gene_index]]).toarray().reshape(-1)
        if toggle:
            # Adding expression data
            plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
            plot_data['Expression'] = expression
            # Subsetting cells
            plot_data = plot_data.loc[cell_filter]
            # Sort values
            plot_data.sort_values(by="Expression", axis = 0, ascending = True, inplace = True)
            # Plotting
            plt.figure(figsize=(6, 5))
            if plot_data.Expression.min() != plot_data.Expression.max():
                ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Expression', hue_norm=(0, plot_data.Expression.max()), palette='viridis', marker='.', s=dot_size, linewidth=0)
                norm = Normalize(0, plot_data.Expression.max())
                sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=ax)
                plt.legend().remove()
                plt.xlabel('UMAP 1')
                plt.ylabel('UMAP 2')
                plt.savefig(f'{self.plot_dir}/{target_gene}_Expression_UMAP.png', bbox_inches='tight', dpi=300)
                plt.close()
            else:
                seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Expression', palette='viridis', marker='.', s=dot_size, linewidth=0)
                legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title=f'{target_gene} Expression', markerscale=markerscale)
                plt.xlabel('UMAP 1')
                plt.ylabel('UMAP 2')
                plt.savefig(f'{self.plot_dir}/{target_gene}_Expression_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
                plt.close()
        else:
            self.log_trace(f'Target gene {target_gene}, not found')
    
    ### ------------------------------------ ###
    
    def plot_gene_set_score(self, set_score: np.ndarray, gene_set_name: str='', datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Check gene set name
        if gene_set_name == '':
            gene_set_name = 'Gene_Set'

        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters) or not hasattr(self, 'clusters'):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Adding expression data
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['Score'] = list(set_score)
        
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
        
        # Sort values
        plot_data.sort_values(by="Score", axis = 0, ascending = True, inplace = True)
        
        # Plotting
        plt.figure(figsize=(6, 5))
        if plot_data.Score.min() != plot_data.Score.max():
            # Set colormap range
            vmin = round(plot_data.Score.min()) if plot_data.Score.min() < 0 else -0.25
            vmax = round(plot_data.Score.max()) if plot_data.Score.max() > 0 else 0.25
            abs_smallest = min(abs(vmin), vmax)
            vmin, vmax = -abs_smallest, abs_smallest
            # Plot
            norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Score', hue_norm=norm, palette='RdBu_r', marker='.', s=dot_size, linewidth=0)
            sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
            plt.legend().remove()
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(f'{self.plot_dir}/{gene_set_name}_Expression_UMAP.png', bbox_inches='tight', dpi=300)
            plt.close()
        else:
            seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Score', palette='viridis', marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Gene Set Expression', markerscale=markerscale)
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(f'{self.plot_dir}/{gene_set_name}_Expression_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
            plt.close()
    
    ### ------------------------------------ ###
    
    def plot_marker_set(self, cell_type_markers: pd.DataFrame, prefix: str='markers', fig_x_size: float=10, fig_y_size: float=6, width_ratios: list[float]=[6, 0.5, 4]) -> None:

        # Format list of markers, their annotation, and their index
        min_markers = 1
        all_markers, markers_annot, all_markers_idx = [], [], []
        for _,(cell_name,markers) in cell_type_markers.iterrows():
            markers_cleaned = [m for m in markers.split(';') if m in self.all_genes]
            markers_idx = [self.all_genes.index(mc) for mc in markers_cleaned]
            if len(markers_cleaned) >= min_markers:
                all_markers += markers_cleaned
                markers_annot += [cell_name for _ in range(len(markers_cleaned))]
                all_markers_idx += markers_idx

        # Format for plotting
        colors = 4 * seaborn.color_palette('tab20')
        unique_markers = np.sort(np.unique(markers_annot)).tolist()
        annot_data = pd.DataFrame({'annotation' : markers_annot,
                                   'color' : [colors[unique_markers.index(ma)] for ma in markers_annot]},
                                  index=all_markers)
        plot_data = pd.DataFrame([],
                                 index=all_markers)
                
        # Add mean expression values
        for cl in np.sort(np.unique(self.clusters)):
            cl_cells = np.where(self.clusters == cl)[0]
            cl_mean = np.mean(csr_array(self.all_data[cl_cells,][:, all_markers_idx]).toarray(), axis=0).tolist()
            plot_data.loc[:, f'cluster_{cl}'] = cl_mean.copy()

        # Remove genes that are not expressed in any sample
        annot_data = annot_data.loc[plot_data.sum(axis=1) > 0]
        plot_data = plot_data.loc[plot_data.sum(axis=1) > 0]

        # Clustering
        cm = seaborn.clustermap(plot_data, z_score=0, center=0, vmin=-3, vmax=3, row_colors=annot_data['color'].values, row_cluster=False, col_cluster=True)
        plt.close()
        col_order = cm.dendrogram_col.reordered_ind # type: ignore
        plot_data = plot_data.iloc[:, col_order]

        # Z score
        plot_data_zscored = zscore(plot_data, axis=1)

        # Color palette
        palette = seaborn.diverging_palette(h_neg=260, h_pos=15, s=75, l=50, sep=1, as_cmap=True)
        # Init plot
        fig, axes = plt.subplots(1, 3, figsize=(fig_x_size, fig_y_size), gridspec_kw={'width_ratios': width_ratios})
        # Heatmap
        seaborn.heatmap(plot_data_zscored, center=0, vmin=-3, vmax=3, cmap=palette, ax=axes[0], cbar_ax=axes[1], cbar=True)
        axes[0].set_yticks(np.arange(0.5, annot_data.shape[0], 1), annot_data.index.values)
        axes[0].tick_params(axis='y', which='major', pad=25, length=0, labelsize=6)
        axes[0].tick_params(axis='x', which='major', labelsize=6)
        # Row annotations
        for i,color in enumerate(annot_data['color'].values):
            axes[0].add_patch(Rectangle(xy=(-0.065, i), width=0.05, height=1, color=color, lw=0,
                                        transform=axes[0].get_yaxis_transform(), clip_on=False))
            axes[0].add_patch(Rectangle(xy=(-0.065, i), width=0.05, height=1, edgecolor='black', facecolor='none', lw=1,
                                        transform=axes[0].get_yaxis_transform(), clip_on=False))

        # Heatmap border
        axes[0].axhline(y=0, color='black', linewidth=2)
        axes[0].axhline(y=plot_data_zscored.shape[0], color='black', linewidth=2)
        axes[0].axvline(x=0, color='black', linewidth=2)
        axes[0].axvline(x=plot_data_zscored.shape[1], color='black', linewidth=2)
        # Patches for the markers annotation
        cell_type_legend_patches = [mpatches.Patch(color=color, label=annot)
                                    for _,(annot,color) in annot_data.drop_duplicates(['annotation', 'color']).iterrows()]
        axes[2].legend(handles=cell_type_legend_patches, title='Marker type', loc='center')
        axes[2].set_axis_off()
        # Save
        plt.tight_layout()
        plt.savefig(f'{self.plot_dir}/{prefix}.png', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_trajectories(self, datasets: list[str]=[], branches: list=[], pseudotime: bool=False, dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()
        
        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        
        # Get list of clusters used in the trajectory
        if not len(branches):
            branches = self.branches
        
        # Get list of clusters used in the trajectory
        clusters = list({vertex for br in branches for vertex in br})
        clusters.sort()
        cell_filter = list(dataset_filter & np.isin(self.clusters, clusters))
        
        # Subsetting cluster color palette
        color_palette = [self.cluster_colors[cl] for cl in clusters]

        # Splitting legend into legend_cols columns
        legend_cols = int(len(clusters) / 16) + 1
        
        # Adding cluster and pseudotime data, then sorting by smallest value
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data['Clusters'] = list(self.clusters)
        if pseudotime:
            plot_data['Pseudotime'] = list(self.pseudotimes)
    
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
    
        # Finding cluster centers for plotting branches
        centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
    
        # Plotting
        if pseudotime:
            # Sort data
            plot_data.sort_values(by="Pseudotime", axis = 0, ascending = True, inplace = True)
            # Plot
            plt.figure(figsize=(6, 5))
            ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Pseudotime', hue_norm=(0, plot_data.Pseudotime.max()), palette='viridis', marker='.', s=dot_size, linewidth=0)
            norm = Normalize(0, plot_data.Pseudotime.max())
            sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
            plt.legend().remove()
            for br in branches:
                cluster_indexes = [clusters.index(cl) for cl in br]
                x, y = centers[cluster_indexes, 0], centers[cluster_indexes, 1]
                plt.plot(x, y, markersize=3, marker='o', linewidth=1, color="red", solid_capstyle='round')
            plt.xlabel('Umap 1')
            plt.ylabel('Umap 2')
            plt.savefig(f'{self.plot_dir}/Trajectories_UMAP_Pseudotime.png', bbox_inches='tight', dpi=300)
        else:
            plt.figure(figsize=(5, 5))
            seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Clusters', palette=color_palette, marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Clusters', ncol=legend_cols, markerscale=markerscale)
            for br in self.branches:
                cluster_indexes = [clusters.index(cl) for cl in br]
                x, y = centers[cluster_indexes, 0], centers[cluster_indexes, 1]
                plt.plot(x, y, markersize=3, marker='o', linewidth=1, color="black", solid_capstyle='round')
            plt.xlabel('Umap 1')
            plt.ylabel('Umap 2')
            plt.savefig(f'{self.plot_dir}/Trajectories_UMAP_NoPseudotime.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
        
    ### ------------------------------------ ###
        
    def plot_variable(self, variable_values: np.ndarray, variable_name: str='Custom variable', categorical: bool=False, sort_values: bool=False, datasets: list[str]=[], clusters: list[int]=[], dot_size: float=1.5, markerscale: float=2) -> None:
        
        # Check plot directory
        self.check_plot_dir()

        # Load PCA embdeddings
        _, _, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Creating a filter for cells of interest
        if not len(datasets):
            dataset_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            dataset_filter = np.array([True if sum([1 for ds in datasets if cell.endswith(f'_{ds}')]) > 0 else False for cell in self.all_cells])
        if not len(clusters) or not hasattr(self, 'clusters'):
            cluster_filter = np.array([True for _ in range(self.all_data.shape[0])])
        else:
            cluster_filter = np.isin(self.clusters, clusters)
        cell_filter = list(dataset_filter & cluster_filter)
        
        # Adding variable data
        plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
        plot_data[variable_name] = list(variable_values)
        
        # Subsetting cells
        plot_data = plot_data.loc[cell_filter]
        
        # Sort values
        if sort_values:
            plot_data.sort_values(by=variable_name, axis=0, ascending=True, inplace=True)
        
        # Plotting
        plot_name = f'{self.plot_dir}/{variable_name.lower().replace(" ", "_")}_UMAP.png'
        plt.figure(figsize=(6, 5))
        if categorical:
            # Subsetting cluster color palette
            color_palette = self.cluster_colors[:len(np.unique(variable_values))]
            seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue=variable_name, palette=color_palette, marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title=variable_name, markerscale=markerscale)
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(plot_name, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
            plt.close()
        elif plot_data[variable_name].min() == plot_data[variable_name].max():
            seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue=variable_name, palette='viridis', marker='.', s=dot_size, linewidth=0)
            legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title=variable_name, markerscale=markerscale)
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(plot_name, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
            plt.close()
        else:
            
            ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue=variable_name, hue_norm=(plot_data[variable_name].min(), plot_data[variable_name].max()), palette='viridis', marker='.', s=dot_size, linewidth=0)
            norm = Normalize(plot_data[variable_name].min(), plot_data[variable_name].max())
            sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
            plt.legend().remove()
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(plot_name, bbox_inches='tight', dpi=300)
            plt.close()

### ---------------------------------------- ###

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

### ---------------------------------------- ###

class ScTypeLike:
    
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
    
    def __init__(self, counts: csr_array, genes: list[str], markers: pd.DataFrame, tissues: list[str]):
        
        # Process markers
        tissue_markers, markers_specificity = self.parse_markers(markers, tissues, genes)
        self.tissue_markers = tissue_markers
        self.markers_specificity = markers_specificity
        
        # Process counts
        weighted_counts = self.parse_counts(counts, genes, markers_specificity)
        self.weighted_counts = weighted_counts
    
    ### ------------------------------------ ###
    ### PARSE MARKERS TABLE                  ###
    ### ------------------------------------ ###
    
    @staticmethod
    def parse_markers(mrks: pd.DataFrame, t: list[str], detected_gns: list[str]) -> tuple[pd.DataFrame, dict[str, float]]:
        
        # Replace nan values
        mrks = mrks.fillna('NONE')
        
        # Subset by tissue
        mrks_sub = mrks.loc[mrks['tissue'].isin(t)]
        
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
    
    def parse_counts(self, cnts: csr_array, gns: list[str], m_specificity: dict[str, float]) -> pd.Series:
        
        # Subset counts matrix
        markers_idx = [gns.index(m) for m in m_specificity.keys()]
        cnts = csr_array(cnts[:, markers_idx])
        
        # Scale counts matrix
        # N.B. will return a full matrix
        scaled_cnts = self.scale_features(cnts)

        # Annotate scaled_cnts
        scaled_cnts = pd.DataFrame(scaled_cnts, columns=list(m_specificity.keys()))
        
        # Weight gene expression
        weighted_cnts = scaled_cnts * pd.Series(m_specificity)
        
        return weighted_cnts

    ### ------------------------------------ ###
    ### ASSIGN CELL TYPES                    ###
    ### ------------------------------------ ###
    
    def score_cell_identity(self) -> None:
        
        # Init score matrix
        scores = pd.DataFrame(np.zeros((self.weighted_counts.shape[0], self.tissue_markers.shape[0])), columns=self.tissue_markers['cell_type'].values)
        
        # Score individual cell types
        for _,(_,cell_type,pos,neg) in self.tissue_markers.iterrows():
            # Score positive markers
            pos = [p for p in pos.split(',') if p in self.markers_specificity.keys()]
            if not len(pos):
                pos_score = np.zeros(scores.shape[0])
            else:
                pos_score = self.weighted_counts[pos].sum(axis=1).values / np.sqrt(len(pos))  # type: ignore
            # Score negative markers
            neg = [n for n in neg.split(',') if n in self.markers_specificity.keys()]
            if not len(neg):
                neg_score = np.zeros(scores.shape[0])
            else:
                neg_score = self.weighted_counts[neg].sum(axis=1).values / np.sqrt(len(neg))  # type: ignore
            # Final score
            cell_type_score = pos_score - neg_score
            scores.loc[:, cell_type] = cell_type_score.copy()
        self.cell_scores = scores
    
        # Score identity based on highest value
        self.cell_identities = scores.idxmax(axis=1).values
    
    ### ------------------------------------ ###
    
    def reduce_to_cluster_level(self, clusters: list[str]) -> None:
        
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
    
    def scale_features(self, mtx: csr_array, gene_mean: np.ndarray=np.array([]), gene_std: np.ndarray=np.array([])) -> np.ndarray:
        
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
    def get_gene_mean_and_std(mtx: csr_array) -> tuple[np.ndarray, np.ndarray]:
        
        # Transposing matrix from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_matrix slicing is efficient along rows
        mtx_t = mtx.transpose()
        
        # Calculate features mean and variance
        # Remember that variance = mean of squared values - square of mean values
        mtx_squared = mtx_t.copy()
        mtx_squared.data **= 2
        gene_mean = np.asarray(mtx_t.mean(axis=1)).reshape(-1)
        gene_std = (np.asarray(mtx_squared.mean(axis=1)).reshape(-1) - (gene_mean ** 2))**0.5
        
        # Reshaping mean and std
        gene_mean = gene_mean.reshape((1, gene_mean.shape[0]))
        gene_std = gene_std.reshape((1, gene_std.shape[0]))
        
        return gene_mean, gene_std
