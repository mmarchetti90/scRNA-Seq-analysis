#!/usr/bin/env python3

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
    10 - generate_psudotime;
    11 - score_cell_cycle or score_gene_set;
"""


### ---------------------------------------- ###

class integrated_analysis:
    
    def __init__(self, max_cells_in_memory=100000, desired_feature_type='Gene Expression'):
        
        # Creating color palette for up to 60 clusters
        """
        # Use the following if there's a lot of clusters
        cluster_colors = [(r/3, g/3, b/3) for r in range(4) for g in range(4) for b in range(4)][1:-1]
        cluster_colors = [cluster_colors[i1 + i2] for i2 in range(4) for i1 in range(0, len(cluster_colors), 4) if i1 + i2 < len(cluster_colors)] # Spacing colors
        random_indexes = [33, 18, 41, 43, 34, 8, 56, 31, 16, 25, 2, 7, 46, 36, 5, 22, 11, 26, 53, 9, 14, 19, 58, 23, 1, 42, 55, 61, 35, 49, 6, 15, 44, 30, 47, 38, 37, 54, 60, 13, 51, 39, 48, 50, 3, 10, 21, 52, 24, 57, 28, 27, 32, 17, 20, 40, 12, 29, 59, 45]
        cluster_colors = [cluster_colors[ri] for ri in random_indexes]
        self.cluster_colors = seaborn.color_palette(cluster_colors, as_cmap=True)
        """
        self.cluster_colors = seaborn.color_palette("tab20")
        
        # Setting maximum number of cells from sparse matrices to load in memory
        self.max_cells_in_memory = max_cells_in_memory
        
        # Setting the selected feature type ('Gene Expression' for RNA data or 'Peak' for ATAC data)
        self.desired_feature_type = desired_feature_type
    
    ### ------------------------------------ ###
    ### DATA LOADING                         ###
    ### ------------------------------------ ###
    
    def load_all_matrices_in_path(self, data_path):
        
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
    
    def load_manifest(self, manifest_path):
        
        try:
            
            self.manifest = pd.read_csv(manifest_path, sep='\t')[['data_label', 'dataset_path', 'normalization']]
            
        except:
            
            print('ERROR: could not read file.')

        if len(np.unique(self.manifest.data_label.values)) < self.manifest.shape[0]:

            print('ERROR: labels are not unique!')
            self.manifest = None
    
    ### ------------------------------------ ###
    
    def load_sparse_matrix(self, main_dir, save_npz_if_missing=True):
        
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
            
                # Constructing a new csr_matrix
                
                matrix_path = matrix_path[0]
                
                # Check matrix files, finding best way to break it into manageable chunks of self.max_cells_in_memory
                matrix_separator, matrix_breaks = self.check_matrix_file(f'{main_dir}/{matrix_path}', cells_num, genes_num, self.max_cells_in_memory)
            
                # Load matrix as a cells-by-genes scipy.sparse.csr_matrix
                sparse_matrix = dok_matrix((cells_num, genes_num))
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

                # Convert to csr_matrix
                sparse_matrix = csr_matrix(sparse_matrix)
                
                # Save as npz file
                if save_npz_if_missing:
                
                    save_npz(f'{main_dir}/{matrix_path}.npz', sparse_matrix, compressed=True)
            
            return barcodes, features, sparse_matrix
            
        except:
            
            return [], [], []
    
    ### ------------------------------------ ###

    def save_sparse_matrix_to_path(self, main_path, label, barcodes, features, matrix, save_mtx=True):
        
        # Saves a sparse matrix in nzip and 10X Genomics format to a desired path after preprocessing
        # N.B. 10X Genomics mtx matrix is saved only if save_mtx = True
        # N.B. The 10X Genomics format was chosen since it's quite human readable and can easily be
        # imported in R Seurat if one desires
        
        # Create dir to store sparse matrix
        out_dir, count = f'{label}', 0
        while out_dir in listdir(main_path):
            
            count += 1
            out_dir = f'{label}_{count}'
        
        makedirs(f'{main_path}/{out_dir}')
        
        # Save barcodes
        out_text = ('\n'.join(barcodes)).encode()
        out_text += b'\n'
        with gzip.open(f'{main_path}/{out_dir}/{label}_barcodes.tsv.gz', 'wb') as output:
            
            output.write(out_text)
        
        # Save features
        out_text = ('\n'.join(features)).encode()
        out_text += b'\n'
        with gzip.open(f'{main_path}/{out_dir}/{label}_features.tsv.gz', 'wb') as output:
            
            output.write(out_text)
        
        # Save matrix as npz
        save_npz(f'{main_path}/{out_dir}/{label}_matrix.mtx.npz', matrix, compressed=True)

        # Save matrix as 10X Genomics mtx file
        if save_mtx:
        
            self.csr_to_mtx(matrix, out_name=f'{main_path}/{out_dir}/{label}_matrix.mtx.gz', max_cells=self.max_cells_in_memory)
        
        return f'{main_path}/{out_dir}'

    ### ------------------------------------ ###
    
    @staticmethod
    def check_matrix_file(path, cells_num, genes_num, max_cells):
        
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
                
                if type(line) == bytes:
                    
                    line = line.decode('utf8')
                
                separator = '\t' if '\t' in line else ' '
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
    def csr_to_mtx(csr, out_name, max_cells=1000):
        
        # Save csr_matrix as 10X Genomics mtx file (100 cells at a time to save memory)
        with gzip.open(out_name, 'wb') as mtx_out:

            # Write header
            out_text = '\n'.join(["%%MatrixMarket matrix coordinate numeric general",
                                  f'{csr.shape[1]}\t{csr.shape[0]}\t{csr.count_nonzero()}'])
            out_text = out_text.encode()
            out_text += b'\n'
            mtx_out.write(out_text)
            
            for i in range(csr.shape[0] // max_cells + 1):
        
                # Creating indexes for mtx rows (genes) and columns (cells)
                # Note that in the csr, rows and columns are transposed
                cols, rows = np.where(csr[i * max_cells : (i + 1) * max_cells].toarray() != 0)
                cols += i * max_cells
                
                # Get values
                counts = csr[cols, rows].A.flatten()
                
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
    def full_to_sparse(full):
        
        # Converts a full matrix to a csr_matrix
        
        # Get features and barcodes
        barcodes = full.columns[1:].to_list()
        features = full['GeneSymbol'].to_list()
        
        # Create counts matrix
        matrix = csr_matrix(full.iloc[:, 1:].T.to_numpy())
        
        return barcodes, features, matrix
    
    ### ------------------------------------ ###
    
    @staticmethod
    def h5_to_csr(h5_path, desired_feature_type=''):

        # Reading HDF5 data
        h5_data = h5py.File(h5_path, 'r')
        
        # Unique cell identifiers
        barcodes = np.array(h5_data['matrix']['barcodes']).astype(str)
        cells_num = len(barcodes)
        
        # Genes
        #gene_ids = np.array(h5_data['matrix']['features']['id']).astype(str)
        gene_names = np.array(h5_data['matrix']['features']['name']).astype(str)
        genes_num = len(gene_names)
        try:
            
            features_type = np.array(h5_data['matrix']['features']['feature_type']).astype(str)
        
        except:
            
            features_type = []
        
        # Init sparse matrix
        sparse_matrix = dok_matrix((cells_num, genes_num))
        
        # Parse sparse
        for N,(start,stop) in enumerate(zip(h5_data['matrix']['indptr'], h5_data['matrix']['indptr'][1:])):

            # Extracting cell data
            counts = h5_data['matrix']['data'][start : stop]
            indexes = h5_data['matrix']['indices'][start : stop]
            
            # Update matrix
            sparse_matrix[N, indexes] = counts
        
        # Filter features based on type
        if len(features_type) and desired_feature_type in ['Gene Expression', 'Peak']:
            
            features_filter = np.where(np.array(features_type) == desired_feature_type)[0]
            
            features = np.array(gene_names)[features_filter].tolist()
        
            sparse_matrix = sparse_matrix[:, features_filter]

        # Convert to csr_matrix
        sparse_matrix = csr_matrix(sparse_matrix)
        
        return barcodes, features, sparse_matrix
    
    ### ------------------------------------ ###
    
    @staticmethod
    def merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=[]):
        
        # Create a metadata file to store dataset labels and number of cells
        metadata = {'datasets_labels' : datasets_labels,
                    'cell_num' : [len(cl) for cl in cells_list]}
        
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
                
                print('ERROR: list not valid')
                return [], [], [], []
        
        else:
            
            # Use genes in common among datasets
            all_genes = np.array([])
            for (index, gl) in enumerate(genes_list):
            
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
            all_data.append(ds[:, genes_index])
        
        all_data = vstack(all_data)
    
        return metadata, all_barcodes, all_genes, all_data
    
    ### ------------------------------------ ###
    
    @staticmethod
    def read_barcodes_file(path):
        
        barcodes = pd.read_csv(path, header=None,  sep='\t').iloc[:, 0].to_list()
        if barcodes[0].lower() in ['barcode', 'barcodes', 'cell', 'cell_id']:
            
            barcodes = barcodes[1:]
        
        return barcodes
    
    ### ------------------------------------ ###
    
    @staticmethod
    def read_features_file(path):
        
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
    def read_matrix_chunk(matrix_path, start, stop, matrix_separator):
        
        # Reads a chunk of a 10X Genomics formatted mtx file
        
        # Define rows to skip and how many to read
        skip_rows = start
        n_rows = stop - start
        
        # Load count matrix
        matrix_chunk = pd.read_csv(matrix_path, header=None, skiprows=skip_rows, nrows=n_rows, skip_blank_lines=True, sep=matrix_separator).to_numpy()
        
        return matrix_chunk
    
    ### ------------------------------------ ###
    ### DATA PREPROCESSING                   ###
    ### ------------------------------------ ###
    
    def normalize_data(self, mtx, skip_first_proportional_fitting=False):
        
        # Normalizes gene counts by proportional fitting + log1p + proportional fitting
        # N.B. first proportional fitting should be skipped if data was already CPM or RPKM normalized
        
        # First proportional fitting
        if not skip_first_proportional_fitting:
            
            mtx = self.proportional_fitting(mtx)
        
        # Log1p transformation
        mtx = mtx.log1p()
        
        # Second proportional fitting
        mtx = self.proportional_fitting(mtx)
        
        return mtx
    
    ### ------------------------------------ ###
    
    def preprocess_samples(self, min_cell_raw_counts=500, min_detected_genes=500):
        
        print('\n########################################')
        print('Preprocessing data')
        
        # Create dir to store sparse matrices after preprocessing
        preprocessed_dir, count = '1_preprocessed_data', 0
        while preprocessed_dir in listdir():
            
            preprocessed_dir = f'1_preprocessed_data_{count}'
            count += 1
        
        mkdir(preprocessed_dir)
        self.preprocessed_dir = preprocessed_dir
        
        # Load and normalize samples, then store them in a sparse matrix
        for (index, row) in self.manifest.iterrows():
            
            label, path, norm_type = row
            
            print('\n########################################')
            print(f'Reading {path}')
            
            # Load data
            if isdir(path): # Data is a sparse matrix
                
                barcodes, features, matrix = self.load_sparse_matrix(path, save_npz_if_missing=False)
            
            elif path.endswith('.h5'): # Data in h5 format
                
                barcodes, features, matrix = self.h5_to_csr(path, self.desired_feature_type)
            
            else: # Data is a full matrix
                
                data = pd.read_csv(path, sep='\t', skip_blank_lines=True)
                print(f'Converting {path} to sparse matrix')
                barcodes, features, matrix = self.full_to_sparse(data)
                
                del data
                gc.collect()
            
            # Skip if dataframe could not be read
            if len(barcodes) == 0:
                
                print(f'WARNING: could not read file {path}')
                continue
            
            # Remove cells with number of genes detected < min_detected_genes
            print('Removing cells with too few genes detected')
            barcodes, matrix = self.remove_cells_with_few_genes(barcodes, matrix, min_detected_genes)
            
            # Normalize data
            if norm_type == 'none': # Performing PFlog1pPF normalization
                
                # Remove cells with read count < min_cell_raw_counts
                print('Removing cells with too few counts')
                barcodes, matrix = self.remove_cells_with_low_counts(barcodes, matrix, min_cell_raw_counts)
                
                # Normalize data
                print('Normalizing counts')
                matrix = self.normalize_data(matrix, skip_first_proportional_fitting=False)
            
            elif norm_type in ['cpm', 'rpkm']: # Performing log1p normalization and proportional fitting
                
                print('Normalizing counts')
                matrix = self.normalize_data(matrix, skip_first_proportional_fitting=True)
            
            elif norm_type in ['pflog1ppf', 'pool_norm']: # Preferred normalizations, data is passed as is
                
                pass
            
            else:
                
                print(f'WARNING: could not determine normalization method of {path}')
                print('Data will be used as is...')
            
            # Append label to cell names
            barcodes = [b + '_' + label for b in barcodes]
            
            # Saving sparse matrix to temporary file
            print('Saving sparse matrix')
            
            _ = self.save_sparse_matrix_to_path(preprocessed_dir, label, barcodes, features, matrix)
            print('########################################\n')
        
        print('Preprocessing done!')
        print(f'Data was saved to sparse matrices in {self.preprocessed_dir}')
        print('########################################\n')

    ### ------------------------------------ ###
    
    @staticmethod
    def proportional_fitting(mtx, target_size=-1):
        
        # Get total cell counts
        total_counts = np.asarray(mtx.sum(axis=1)).reshape(-1)
        
        # Compute correction factors
        if target_size == -1:
        
            target_size = np.median(total_counts)
        
        else:
            
            pass
        
        correction_factors = target_size / total_counts
        correction_factors = correction_factors.reshape((len(correction_factors), 1))
        
        # Correct matrix
        mtx = csr_matrix(mtx.multiply(correction_factors))

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
    
    def find_common_highly_variable_genes(self, data_path, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000, union_mode=False):
        
        # Wrapper for finding highly variable genes (HVGs) common to all datasets
        # If union_mode = True, all HVGs from individual datasets are combined
        
        print('\n########################################')
        print(f'Finding {"all" if union_mode else "common"} highly variable genes')
        
        # Create dir to store HVGs
        hvgs_dir, count = '2_hvgs', 0
        while hvgs_dir in listdir():
            
            hvgs_dir = f'2_hvgs_{count}'
            count += 1
        
        mkdir(hvgs_dir)
        
        # Load matrices
        datasets_labels, datasets, _, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Find variable features in each dataset, then intersect
        hvgs = np.array([])
        
        for label, ds, gl in zip(datasets_labels, datasets, genes_list):
            
            print('\n########################################')
            print(f'Finding highly variable genes for {label}')
            
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

            print(f'Found {len(dataset_hvgs)} highly variable genes')            
            print('########################################\n')
        
        # Store HVGs
        self.hvgs = list(hvgs)
        
        # Save HVGs to file
        out_text = ('\n'.join(self.hvgs)).encode()
        out_text += b'\n'
        with gzip.open(f'{hvgs_dir}/common_higly_variable_genes.tsv.gz', 'wb') as output:
            
            output.write(out_text)
            
        print(f'Found {len(hvgs)} common highly variable genes')
        print(f'Highly variable genes were saved to {hvgs_dir}/common_higly_variable_genes.tsv.gz')
    
    ### ------------------------------------ ###
    
    def find_highly_variable_genes_in_whole_dataset(self, data_path, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000):
        
        # Find highly variable genes after merging all datasets into one
        # Use this if find_common_highly_variable_genes returns very few genes
        
        print('\n########################################')
        print('Finding highly variable genes in merged datasets')
        
        # Create dir to store HVGs
        hvgs_dir, count = '2_hvgs', 0
        while hvgs_dir in listdir():
            
            hvgs_dir = f'2_hvgs_{count}'
            count += 1
        
        mkdir(hvgs_dir)
        
        # Load matrices
        datasets_labels, datasets, _, genes_list = self.load_all_matrices_in_path(data_path)
        
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
        with gzip.open(f'{hvgs_dir}/higly_variable_genes_whole_dataset.tsv.gz', 'wb') as output:
            
            output.write(out_text)
            
        print(f'Found {len(hvgs)} highly variable genes')
        print(f'Highly variable genes were saved to {hvgs_dir}/higly_variable_genes_whole_dataset.tsv.gz')
        print('########################################\n')
    
    ### ------------------------------------ ###
    
    def load_hvgs(self, data_path):
        
        try:
            
            self.hvgs = pd.read_csv(data_path, sep='\t', header=None).to_numpy().flatten().tolist()
        
        except:
            
            print("ERROR: couldn't read file")
            self.hvgs = []
    
    ### ------------------------------------ ###
    
    @staticmethod
    def find_highly_variable_genes(mtx, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000):
        
        # Graph based, similar to Seurat's FindVariableGenes
        
        # Transposing matrix from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_matrix slicing is efficient along rows
        mtx = mtx.transpose()
        
        # Calculate features mean and variance
        # Remember that variance = mean of squared values - square of mean values
        print('Computing genes mean and variance')
        mtx_squared = mtx.copy()
        mtx_squared.data **= 2
        gene_mean = np.asarray(mtx.mean(axis=1)).reshape(-1)
        gene_variance = np.asarray(mtx_squared.mean(axis=1)).reshape(-1) - (gene_mean ** 2)
        
        # Calculating e^mean and log(variance/mean) (X and Y respectively) for each feature
        var_plot = pd.DataFrame({'GeneID' : range(0, mtx.shape[0]),
                                 'X' : np.exp(gene_mean),
                                 'Y' : np.log(gene_variance / (gene_mean + 0.00001))})
        min_x, max_x = max(x_low_cutoff, var_plot.X.min()), min(x_high_cutof, var_plot.X.max())
        bin_width = (max_x - min_x) / nbins

        # Binning and computing log(variance/mean) z-scores
        print('Binning genes and computing log(variance/mean) z-scores')
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
        most_variable_features = list(most_variable_features.loc[(most_variable_features.ZScores >= y_low_cutoff) &
                                                                 (most_variable_features.ZScores <= y_high_cutoff),].iloc[:max_features,].GeneID)
        
        return most_variable_features
    
    ### ------------------------------------ ###
    ### DATA INTEGRATION                     ###
    ### ------------------------------------ ###
    
    def integrate_data(self, data_path='', use_highly_variable_genes=True, union=False, adjust_library_size=True, use_scanorama_embeddings=True, max_pca_components=50):
        
        # Data integration with Scanorama
        
        print('\n########################################')
        print('Integrating data using Scanorama')
        
        # Define path to preprocessed sparce matrices folder
        if not len(data_path):
            
            try:
                
                data_path = self.preprocessed_dir
                
            except:
                
                print('ERROR: missing data_path value')
                return ''
        
        # Check that data_path exists
        if not exists(data_path):
            
            print('ERROR: bad data_path value')
            return ''
        
        # Create directory to store the integrated sparse matrix
        integrated_data_dir, count = '3_integrated_data', 0
        while integrated_data_dir in listdir():
            
            integrated_data_dir = f'3_integrated_data_{count}'
            count += 1
        
        mkdir(integrated_data_dir)
        self.integrated_data_dir = integrated_data_dir
        
        # Load matrices for Scanorama
        datasets_labels, datasets, cells_list, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Correcting for library size
        if adjust_library_size:
        
            datasets_target_size = np.median(np.concatenate([ds.sum(axis=1).A.ravel() for ds in datasets]))
            datasets = [self.proportional_fitting(ds, datasets_target_size) for ds in datasets]
        
        # Filter for common HVGs
        if use_highly_variable_genes:
            
            # Merging datasets using HVGs, then splitting into individual datasets. This so as to filter datasets for HVGs
            datasets_metadata, _, all_genes, all_data = self.merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=self.hvgs)
            datasets_metadata = pd.DataFrame(datasets_metadata)
            
            # Splitting all_data into its constituent datasets
            for index in range(datasets_metadata.shape[0]):
                
                # Get dataset info
                label = datasets_metadata.iloc[index, 0]
                start_row = datasets_metadata.iloc[:index, 1].values.sum()
                end_row = datasets_metadata.iloc[:index + 1, 1].values.sum()
                
                datasets[index] = all_data[start_row : end_row]
            
            genes_list = [all_genes for _ in range(len(datasets))]
            
        # Integrate data
        if not use_scanorama_embeddings:
            
            corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=False, ds_names=datasets_labels, approx=True, knn=20, hvg=None, return_dense=False, union=union, seed=42, batch_size=self.max_cells_in_memory)
            
        else:
            
            integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True, ds_names=datasets_labels, approx=True, knn=20, hvg=None, return_dense=False, union=union, seed=42, batch_size=self.max_cells_in_memory)
            
            integrated = np.concatenate(integrated, axis=0)[:, :max_pca_components]
            
            self.pca_data = integrated
            
            pd.DataFrame(integrated).to_csv(f'{integrated_data_dir}/integrated_pca_embdeddings.tsv.gz', sep='\t', header=False, index=False)
        
        # Save matrices
        for (label, matrix, barcodes, features) in zip(datasets_labels, datasets, cells_list, genes_list):
            
            print(f'Saving {label} corrected matrix to file')
            self.save_sparse_matrix_to_path(integrated_data_dir, label, barcodes, features, matrix, save_mtx=True)
        
        print(f'All done!\nData was saved to sparse matrices in {self.integrated_data_dir}')
        print('########################################\n')
        
        # Merging integrated dataset
        self.merge_integrated_datasets(self.integrated_data_dir, False)
    
    ### ------------------------------------ ###
    
    def merge_integrated_datasets(self, data_path='', adjust_library_size=False):
        
        print('\n########################################')
        print('Merging integrated data to single dataset')
        
        # Define path to integrated sparce matrices folder
        if not len(data_path):
            
            try:
                
                data_path = self.integrated_data_dir
                
            except:
                
                print('ERROR: missing data_path value')
                return ''
        
        # Check that data_path exists
        if not exists(data_path):
            
            print('ERROR: bad data_path value')
            return ''
        
        # Create directory to store the merged sparse matrix
        merged_data_dir, count = '4_merged_data', 0
        while merged_data_dir in listdir():
            
            merged_data_dir = f'4_merged_data_{count}'
            count += 1
        
        mkdir(merged_data_dir)
        self.merged_data_dir = merged_data_dir
        
        # Load integrated data
        datasets_labels, datasets, cells_list, genes_list = self.load_all_matrices_in_path(data_path)
        
        # Correcting for library size
        if adjust_library_size:
        
            datasets_target_size = np.median(np.concatenate([ds.sum(axis=1).A.ravel() for ds in datasets]))
            datasets = [self.proportional_fitting(ds, datasets_target_size) for ds in datasets]
        
        # Merge datasets
        datasets_metadata, all_cells, all_genes, all_data = self.merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=[])
        datasets_metadata = pd.DataFrame(datasets_metadata)
        
        # Store in class
        self.datasets_metadata = datasets_metadata
        self.all_cells = all_cells
        self.all_genes = all_genes
        self.all_data = all_data
        
        # Save metadata
        print('Saving datasets metadata to file')
        datasets_metadata.to_csv(f'{merged_data_dir}/datasets_metadata.tsv.gz', sep='\t', index=None)

        # Save barcodes
        print('Saving merged barcodes to file')
        out_text = ('\n'.join(all_cells)).encode()
        out_text += b'\n'
        with gzip.open(f'{merged_data_dir}/all_barcodes.tsv.gz', 'wb') as output:
            
            output.write(out_text)
        
        # Save features
        print('Saving merged gene names to file')
        out_text = ('\n'.join(all_genes)).encode()
        out_text += b'\n'
        with gzip.open(f'{merged_data_dir}/all_features.tsv.gz', 'wb') as output:
            
            output.write(out_text)
        
        # Save matrix as npz
        print('Saving integrated sparse matrix to file')
        save_npz(f'{merged_data_dir}/all_data_matrix.mtx.npz', all_data, compressed=True)

        # Save matrix as 10X Genomics mtx file
        print('Saving integrated sparse matrix to file in 10X Genomics format')
        self.csr_to_mtx(all_data, out_name=f'{merged_data_dir}/all_data_matrix.mtx.gz', max_cells=self.max_cells_in_memory)
        
        print(f'All done!\nData was saved to sparse matrices in {self.merged_data_dir}')
        print('########################################\n')
    
    ### ------------------------------------ ###
    
    def load_merged_datasets(self, data_path):
        
        try:
            
            self.datasets_metadata = pd.read_csv(f'{data_path}/datasets_metadata.tsv.gz', sep='\t', header=0)
            self.all_cells, self.all_genes, self.all_data = self.load_sparse_matrix(data_path, save_npz_if_missing=False)
            self.merged_data_dir = data_path
            
        except:
            
            print("ERROR: couldn't read files")
            self.merged_data_dir = ''
    
    ### ------------------------------------ ###
    ### DIMENSIONALITY REDUCTION             ###
    ### ------------------------------------ ###
    
    def reduce_dimensions(self, use_integrated_pca=True, integrated_pca_path='', use_highly_variable_genes=True, max_pca_components=50, neighbors=30, n_permutations=3):
        
        # Wrapper for umap_embedding and pca_umap_embdedding functions
        # N.B. PCA/UMAP models and embeddings are not kept in memory, but have to be loaded with self.load_reduced_dimensions afterwards
        
        if use_integrated_pca:
            
            self.umap_embedding(integrated_pca_path, neighbors, n_permutations)
        
        else:
            
            self.pca_umap_embdedding(use_highly_variable_genes, max_pca_components, neighbors, n_permutations)
    
    ### ------------------------------------ ###
    
    def pca_umap_embdedding(self, use_highly_variable_genes, max_pca_components, neighbors, n_permutations):
        
        print('\n########################################')
        print('Reducing dataset dimensions')
        
        # Create dir to store dimensionality reduction data
        lower_dimensions_dir, count = '5_lower_dimensions', 0
        while lower_dimensions_dir in listdir():
            
            lower_dimensions_dir = f'5_lower_dimensions_{count}'
            count += 1
        
        mkdir(lower_dimensions_dir)
        self.lower_dimensions_dir = lower_dimensions_dir
        
        # Filter for HVGs
        all_data = self.all_data.copy()
        if use_highly_variable_genes:
            
            all_data = all_data[:, [self.all_genes.index(hvg) for hvg in self.hvgs]]
            
        # Create subsets of self.max_cells_in_memory cells to fit the PCA and UMAP models to
        # N.B. All cells are passed through 3 times in different groupings
        if all_data.shape[0] < self.max_cells_in_memory:
            
            cell_subsets = self.make_cell_subsets(all_data.shape[0], self.max_cells_in_memory)
        
        else:
            
            cell_subsets = [subset for _ in range(n_permutations) for subset in self.make_cell_subsets(all_data.shape[0], self.max_cells_in_memory)]
            
        # Calculate gene means and std across all cells to be used for scaling
        gene_mean, gene_std = self.get_gene_mean_and_std(all_data)
        
        # Fit PCA embedding using randomized subsets of all_data
        print('Fitting PCA model')
        pca_model = PCA(max_pca_components)
        for (i, cs) in enumerate(cell_subsets):
            
            # Scale features for a subset of cells and fit PCA
            all_data_subset = self.scale_features(all_data[cs,], gene_mean, gene_std)
            
            # Fit PCA model
            pca_model.fit(all_data_subset)
    
        # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
        x0, x1 = 0, min(len(pca_model.explained_variance_), max_pca_components - 1)
        y0, y1 = pca_model.explained_variance_[x0], pca_model.explained_variance_[x1]
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(pca_model.explained_variance_[:max_pca_components])]
        optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
        print(f'Optimal PCA components: {optimal_pca_components}')
    
        # Elbow plot of explained variance
        plt.figure(figsize=(5, 3))
        plt.plot(range(x0 + 1, x1 + 2), pca_model.explained_variance_[:max_pca_components] / 100, 'b', marker='o', markersize=5, linewidth=1)
        plt.plot([optimal_pca_components, optimal_pca_components], [y0 / 100, y1 / 100], linestyle='dashed', color='red', linewidth=1)
        plt.xlabel('PC')
        plt.ylabel('Explained Variance (%)')
        plt.tight_layout()
        plt.savefig(f'{lower_dimensions_dir}/PCA_ExplainedVariance', dpi=300)
        plt.close()
        
        # Fit UMAP embedding using randomized subsets of all_data
        print('Fitting UMAP model')
        umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
        for (i, cs) in enumerate(cell_subsets):
        
            # Scale features for a subset of cells
            all_data_subset = self.scale_features(all_data[cs,], gene_mean, gene_std)
            
            # PCA transform
            all_data_subset = pca_model.transform(all_data_subset)
        
            # Fit UMAP model
            umap_model.fit(all_data_subset[:, :optimal_pca_components])
        
        # Save gene means and std used for scaling (useful to add new cells if one does not want to refir PCA and UMAP)
        pd.DataFrame(gene_mean).to_csv(f'{lower_dimensions_dir}/gene_means.tsv.gz', sep='\t', header=False, index=False)
        pd.DataFrame(gene_std).to_csv(f'{lower_dimensions_dir}/gene_stds.tsv.gz', sep='\t', header=False, index=False)
        
        # Save PCA and UMAP models
        print('Saving PCA and UMAP models to pickle file')
        pk.dump(pca_model, open(f'{lower_dimensions_dir}/pca_optimalcomp-{optimal_pca_components}.pkl', 'wb'))
        pk.dump(umap_model, open(f'{lower_dimensions_dir}/umap_neighbors-{neighbors}.pkl', 'wb'))
        
        # Save transformed data
        # N.B Scaled data is saved as a cell-to-gene full matrix, so it will be large
        print('Saving embedded data')
        #scaled_out = gzip.open(f'{lower_dimensions_dir}/scaled_data.tsv.gz', 'wb')
        pca_out = gzip.open(f'{lower_dimensions_dir}/pca_data.tsv.gz', 'wb')
        umap_out = gzip.open(f'{lower_dimensions_dir}/umap_data.tsv.gz', 'wb')
        
        for i in range(all_data.shape[0] // self.max_cells_in_memory + 1):
            
            # Scale features for a subset of cells
            all_data_subset = self.scale_features(all_data[i * self.max_cells_in_memory : (i + 1) * self.max_cells_in_memory,], gene_mean, gene_std)
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
        
        print('All done!')
        print(f'PCA and UMAP models were stored in {lower_dimensions_dir}')
        print('########################################\n')
        
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
    
    def umap_embedding(self, integrated_pca_path, neighbors, n_permutations):
        
        # N.B. PCA/UMAP models and embeddings are not kept in memory, but have to be loaded with self.load_reduced_dimensions afterwards
        
        print('\n########################################')
        print('Reducing dataset dimensions')
        
        # Create dir to store dimensionality reduction data
        lower_dimensions_dir, count = '5_lower_dimensions', 0
        while lower_dimensions_dir in listdir():
            
            lower_dimensions_dir = f'5_lower_dimensions_{count}'
            count += 1
        
        mkdir(lower_dimensions_dir)
        self.lower_dimensions_dir = lower_dimensions_dir
        
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
        print('Fitting UMAP model')
        umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
        for (i, cs) in enumerate(cell_subsets):
        
            # Subset cells
            all_data_subset = pca_data[cs,]
        
            # Fit UMAP model
            umap_model.fit(all_data_subset)
        
        # Save UMAP model
        print('Saving UMAP model to pickle file')
        pk.dump(umap_model, open(f'{lower_dimensions_dir}/umap_neighbors-{neighbors}.pkl', 'wb'))
        
        # Save transformed data
        pd.DataFrame(pca_data).to_csv(f'{lower_dimensions_dir}/pca_data.tsv.gz', sep='\t', header=False, index=False)
        
        print('Saving embedded data')
        umap_out = gzip.open(f'{lower_dimensions_dir}/umap_data.tsv.gz', 'wb')
        
        for i in range(pca_data.shape[0] // self.max_cells_in_memory + 1):
            
            # UMAP transform
            umap_data_subset = umap_model.transform(pca_data[i * self.max_cells_in_memory : (i + 1) * self.max_cells_in_memory,])
            out_text = '\n'.join(['\t'.join(line.astype(str)) for line in umap_data_subset])
            out_text = out_text.encode() + (b'\n' if i < self.all_data.shape[0] // self.max_cells_in_memory else b'')
            umap_out.write(out_text)
        
        # Close files
        umap_out.close()
        
        print('All done!')
        print(f'PCA and UMAP models were stored in {lower_dimensions_dir}')
        print('########################################\n')
    
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
    
    @staticmethod
    def load_reduced_dimensions(lower_dimensions_dir):
        
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
            
            print('ERROR: make sure all the needed files exist in folder.')
            return PCA(), np.array([]), 0, umap.UMAP(), np.array([])
    
    ### ------------------------------------ ###

    @staticmethod
    def make_cell_subsets(N, max_cells):
        
        subsets, indexes = [], np.arange(N)
        
        while len(indexes):
            
            new_subset = np.random.choice(indexes, min(max_cells, len(indexes)), replace=False)
            subsets.append(new_subset)
            indexes = indexes[~ np.isin(indexes, new_subset)]
        
        return subsets

    ### ------------------------------------ ###
    ### CLUSTERING                           ###
    ### ------------------------------------ ###

    def check_cluster_composition(self):
        
        try:
            
            _ = self.clusters
            
        except:
            
            print('ERROR: need to run/load clustering analysis first')
        
        # Get unique clusters
        clusters_id = np.sort(np.unique(self.clusters))
        
        # Get unique datasets
        datasets_labels = self.datasets_metadata.datasets_labels.to_list()
        
        # Get total number of cells per dataset
        datasets_cells = self.datasets_metadata.cell_num.to_list()
        
        # Init dict for results
        results = {key : [] for key in ['cluster', 'cluster_cells'] + datasets_labels + ['hypergeom_test_flag']}
        results['cluster'] = clusters_id
        
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
        
    def cluster_cells(self, n_neighbors=10, resolution_parameter=1.0, beta=0.01):
        
        print('\n########################################')
        print('Clustering cells')
        
        # Create dir to store cluster embeddings
        clusters_dir, count = '6_clustering', 0
        while clusters_dir in listdir():
            
            clusters_dir = f'6_clustering_{count}'
            count += 1
        
        mkdir(clusters_dir)
        self.clusters_dir = clusters_dir
        
        # Setting random seed (helps with clustering consistency)
        random.seed(42)
        
        # Check that self.lower_dimensions_dir exists
        try:
            
            _ = self.lower_dimensions_dir
        
        except:
            
            print('ERROR: missing path to lower dimension data (self.lower_dimensions_dir)')
            return np.array([])
        
        # Load PCA and UMAP embdeddings
        _, pca_data, _, _, umap_data = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Computing kneighbors sparse matrix
        kneighbors_matrix = kneighbors_graph(X=pca_data, n_neighbors=n_neighbors)
        
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
        clusters = graph.community_leiden(objective_function='modularity', weights=None, resolution_parameter=resolution_parameter, beta=beta, initial_membership=None, n_iterations=2, node_weights=None).membership
        clusters = np.array(clusters)
        
        self.clusters = clusters
        
        # Save clusters to file
        pd.DataFrame(clusters).to_csv(f'{clusters_dir}/cluster_embeddings.tsv.gz', sep='\t', header=None, index=None)
        
        # Check cluster composition
        self.check_cluster_composition()
        
        print('All done!')
        print(f'Clusters info was stored in {clusters_dir}/cluster_embeddings.tsv.gz')
        print('########################################\n')
     
    ### ------------------------------------ ###
    
    def find_clusters_markers(self, min_fc=0.25, min_pct=0.1, min_cells_in_cluster=50):
        
        print('\n########################################')
        
        # Create dir to store cluster markers
        cluster_markers_dir, count = '7_cluster_markers', 0
        while cluster_markers_dir in listdir():
            
            cluster_markers_dir = f'7_cluster_markers_{count}'
            count += 1
        
        mkdir(cluster_markers_dir)
        self.cluster_markers_dir = cluster_markers_dir
        
        # Init results table
        cluster_markers = {"cluster" : [], "gene_name" : [], "cluster_pct" : [], "other_clusters_pct" : [], "cluster_expression" : [], "other_clusters_expression" : [], "log2FC" : [], "pval" : [], "padj" : []}
        
        # Get clusters ids
        clusters_id = np.sort(np.unique(self.clusters))
        
        # Find markers in each cluster vs other clusters
        for cl in clusters_id:
            
            cluster_cells = np.where(self.clusters == cl)[0]
            other_cells = np.where(self.clusters != cl)[0]
            
            if len(cluster_cells) < min_cells_in_cluster: # Skip clusters with too few cells
                
                print(f'WARNING: cluster {cl} comprises only {len(cluster_cells)} cells. Skipping cluster...')
                continue
            
            print(f'Finding markers of cluster {cl}')
            
            significant_genes_name, log2_fc, pop1_expression, pop2_expression, pct_pop1, pct_pop2, pvals, padjs = self.compare_populations(self.all_genes, self.all_data, cluster_cells, other_cells, self.max_cells_in_memory, min_fc, min_pct)
            
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
        cluster_markers.to_csv(f'{cluster_markers_dir}/cluster_markers.tsv.gz', sep='\t', header=True, index=False)
        
        print('All done!')
        print(f'Clusters markers were saved to {cluster_markers_dir}/cluster_markers.tsv.gz')
        print('########################################\n')
    
    ### ------------------------------------ ###
    
    def load_clusters(self, data_path):
        
        try:
            
            self.clusters = pd.read_csv(data_path, sep='\t', header=None).to_numpy().flatten()
            self.clusters_dir = '/'.join(data_path.split('/')[:-1])
        
        except:
            
            print("ERROR: couldn't read file")
            self.clusters = np.array([])
            self.clusters_dir = ''
    
    ### ------------------------------------ ###
    
    def load_cluster_identities(self, data_path):
        
        try:
            
            self.cluster_identities = pd.read_csv(data_path, sep='\t', header=0)
        
        except:
            
            print("ERROR: couldn't read file")
            self.cluster_identities = pd.DataFrame()
    
    ### ------------------------------------ ###
    
    def load_cluster_markers(self, data_path):
        
        try:
            
            self.cluster_markers = pd.read_csv(data_path, sep='\t', header=0)
        
        except:
            
            print("ERROR: couldn't read file")
            self.cluster_markers = pd.DataFrame()
    
    ### ------------------------------------ ###
    
    def score_cluster_identity(self, ref_expression, ref_samples, clusters_id=[]):
        
        # Similar to CIPR, but simpler and better
        
        print('\n########################################')
        print('Assigning identities to clusters')
        
        # Create dir to store dimensionality reduction data
        cluster_identity_dir, count = '8_cluster_identities', 0
        while cluster_identity_dir in listdir():
            
            cluster_identity_dir = f'8_cluster_identities_{count}'
            count += 1
        
        mkdir(cluster_identity_dir)
        self.cluster_identity_dir = cluster_identity_dir
        
        # Get clusters ids
        if not len(clusters_id):
            
            clusters_id = np.sort(np.unique(self.clusters))

        # Listing possible cell types and subtypes
        cell_classes = ref_samples.reference_cell_type
        cell_subclasses = ref_samples.long_name
        
        # Init cluster identity matrix
        cluster_identities = {'cluster' : [],
                              'cell_type' : [],
                              'cell_type_prob' : [],
                              'cell_sub_type' : [],
                              'cell_sub_type_prob' : []}
        
        # Assessing identity probabilities for each cluster
        for cl in clusters_id:
            
            print(f'Scoring cluster {cl} identity')
            
            # Subset markers
            cl_markers = self.cluster_markers.loc[(self.cluster_markers.cluster == cl) &
                                                  (self.cluster_markers.padj < 0.05) &
                                                  (self.cluster_markers.gene_name.str.lower().isin(ref_expression.iloc[:, 0].str.lower())),]
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
                
                class_prob = probs[ref_samples.reference_cell_type == c_class].sum()
                
                if class_prob > best_class_prob:
                    
                    best_class = c_class
                    best_class_prob = class_prob

            # Finding most probable subclass
            subclass_probs = probs[ref_samples.reference_cell_type == best_class]
            best_subclass_prob = subclass_probs.max()
            best_subclass = list(cell_subclasses[ref_samples.reference_cell_type == best_class][subclass_probs == best_subclass_prob])[0]
            
            # Adding cluster info
            cluster_identities['cluster'].append(cl)
            cluster_identities['cell_type'].append(best_class)
            cluster_identities['cell_type_prob'].append(best_class_prob)
            cluster_identities['cell_sub_type'].append(best_subclass)
            cluster_identities['cell_sub_type_prob'].append(best_subclass_prob)
        
        cluster_identities = pd.DataFrame(cluster_identities)
        
        self.cluster_identities = cluster_identities
        
        # Saving to file
        cluster_identities.to_csv(f'{cluster_identity_dir}/cluster_identities.tsv.gz', sep='\t', header=True, index=False)
        
        print('All done!')
        print(f'Clusters identities were saved to {cluster_identity_dir}/cluster_identities.tsv.gz')
        print('########################################\n')

    ### ------------------------------------ ###
    
    @staticmethod
    def compare_populations(gene_names, data, pop1, pop2, max_cells_in_memory=10000, min_fc=0.25, min_pct=0.1):
        
        # Split all_data into pop1 and pop2
        pop1_data = data[pop1, ]
        pop2_data = data[pop2, ]
        
        # Transposing matrices from cell-to-genes to genes-to-cells to make calculations (slightly) faster
        # Remember that csr_matrix slicing is efficient along rows
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
            
            pvals = mannwhitneyu(pop1_data.A, pop2_data.A, alternative="two-sided", axis=0).pvalue
            
        else:
            
            pvals = []
            for gene in range(pop1_data.shape[1]):
                
                pval = mannwhitneyu(pop1_data[:, gene].A, pop2_data[:, gene].A, alternative="two-sided", axis=0).pvalue[0]
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
    
    def create_mst(self, clusters_id=[]):
        
        # Kruskal Minimal Spanning Tree algorithm
        
        try:
            
            _ = self.lower_dimensions_dir
            
        except:
            
            print('ERROR: missing path to lower dimension data (self.lower_dimensions_dir)')
        
        print('\n########################################')
        print('Starting trajectory analysis')
        
        # Create dir to store dimensionality reduction data
        trajectories_dir, count = '9_trajectories', 0
        while trajectories_dir in listdir():
            
            trajectories_dir = f'9_trajectories_{count}'
            count += 1
        
        mkdir(trajectories_dir)
        self.trajectories_dir = trajectories_dir
        
        # Load PCA embdeddings
        _, pca_data, _, _, _ = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Finding clusters centroids (graph vertices)
        if not len(clusters_id):
            
            clusters_id = np.sort(np.unique(self.clusters))
        
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
        file_name = self.save_trajectory_branches(trajectories_dir, branches)
        self.trajectory_file_name = file_name
        
        print('All done!')
        print(f'Trajectory analysis was saved to {trajectories_dir}/{file_name}')
        print('########################################\n')

    ### ------------------------------------ ###
    
    def generate_pseudotime(self):
        
        try:
            
            _ = self.lower_dimensions_dir
            
        except:
            
            print('ERROR: missing path to lower dimension data (self.lower_dimensions_dir)')
        
        try:
            
            _ = self.trajectories_dir
            _ = self.trajectory_file_name
            
        except:
            
            print('ERROR: please run or load trajectory analysis')
        
        print('\n########################################')
        print('Calculating cells pseudotime along trajectories')
        
        # Load PCA embdeddings
        _, pca_data, _, _, _ = self.load_reduced_dimensions(self.lower_dimensions_dir)
        
        # Listing clusters used
        clusters_id = list({vertex for br in self.branches for vertex in br})
        clusters_id.sort()
        
        # Calculating cluster centroids
        centroids = np.array([pca_data[self.clusters == cl, ].mean(axis = 0) for cl in clusters_id])
        
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
                cell_pseudotime = ((projection**2).sum())**0.5
                
            elif vertex_index == len(branch) - 1: # Cell belongs to end vertex. Pseudotime = projection on last edge of branch
                
                # Projecting cell to last edge of branch
                branch_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[-2:]])
                edge_vector = np.diff(branch_centroids, axis=0)
                projection = (edge_vector * coords) / ((edge_vector**2).sum())**0.5
          
                # Calculating previous branch length
                previous_centroids = np.array([centroids[clusters_id.index(cl)] for cl in branch[:vertex_index]])
                added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
          
                # Computing pseudotime
                cell_pseudotime = ((projection**2).sum())**0.5 + added_distance
                
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
                added_distance = sum(((np.diff(previous_centroids)**2).sum(axis=1))**0.5)
          
                # Computing pseudotime
                cell_pseudotime = ((projection**2).sum())**0.5 + added_distance
                
            pseudotimes.append(cell_pseudotime)
        
        # Normalize pseudotime to 0-100
        pseudotimes = (np.array(pseudotimes) - min(pseudotimes)) * 100 / (max(pseudotimes) - min(pseudotimes))
        self.pseudotimes = pseudotimes
      
        # Save pseudotime to path
        out_path = f'{self.trajectories_dir}/{self.trajectory_file_name}'.replace('.txt', '_pseudotime.tsv.gz')
        pd.DataFrame(pseudotimes).to_csv(out_path, sep='\t', header=False, index=False)
        
        print('All done!')
        print(f'Trajectory analysis was saved to {out_path}')
        print('########################################\n')
        
    ### ------------------------------------ ###
    
    def load_trajectory_branches(self, data_path):
        
        try:
            
            self.trajectories_dir = '/'.join(abspath(data_path).split('/')[:-1])
            self.trajectory_file_name = data_path.split('/')[-1]
            self.branches = [[int(cluster) for cluster in line.split('\t')] for line in open(data_path).read().split('\n')]
            
        except:
            
            print("ERROR: couldn't read file")
            self.trajectories_dir = '/'.join(abspath(data_path).split('/')[:-1])
            self.trajectory_file_name = ''
            self.branches = []
    
    ### ------------------------------------ ###
    
    def load_trajectory_pseudotime(self, data_path):
        
        try:
            
            self.pseudotimes = pd.read_csv(data_path, sep='\t', header=None).iloc[:, 0].to_numpy()
            
        except:
            
            print("ERROR: couldn't read file")
            self.pseudotimes = np.array([])
    
    ### ------------------------------------ ###
    
    def save_trajectory_branches(self, save_dir, branches):
        
        # Create label to store for trajectories file
        file_name, count = 'trajectory_analysis.txt', 0
        while f'{save_dir}/{file_name}' in listdir(save_dir):
            
            file_name = f'trajectory_analysis_{count}.txt'
            count += 1
        
        # Format output text
        out_text = '\n'.join(['\t'.join([str(b) for b in br]) for br in branches])
        
        # Save trajectory to file
        with open(f'{save_dir}/{file_name}', 'w') as output:
            
            output.write(out_text)
        
        return file_name
    
    ### ------------------------------------ ###
    ### CELL CYCLE ANALYSIS                  ###
    ### ------------------------------------ ###
    
    def score_cell_cycle(self, g1pm_features, s_features, g2m_features, gene_pool=[], bins=50, ctrl_genes_num=50):
        
        print('\n########################################')
        print('Scoring cell cycle')
        
        # Create dir to store dimensionality reduction data
        cell_cycle_dir, count = '10_cell_cycle', 0
        while cell_cycle_dir in listdir():
            
            cell_cycle_dir = f'10_cell_cycle_{count}'
            count += 1
        
        mkdir(cell_cycle_dir)
        self.cell_cycle_dir = cell_cycle_dir
        
        phases = ["G1PM", "S", "G2M"]
        
        # Scoring gene sets
        toggle = 1
        cell_cycle_scores = []
        for set_name, gene_set in zip(phases, [g1pm_features, s_features, g2m_features]):
            
            print(f'Scoring {set_name} phase')
            
            set_score = self.score_gene_set(gene_set, gene_pool, bins, ctrl_genes_num)

            if set_score is None:
                
                toggle = 0
                bad_gene_set = set_name
                break
            
            else:
                
                cell_cycle_scores.append(set_score)
        
        if toggle == 0:
            
            print(f'Error: {bad_gene_set} list is bad')
            return None
        
        else:
            
            cell_cycle_scores = np.stack(cell_cycle_scores).T
            cell_cycle_scores = pd.DataFrame(cell_cycle_scores, columns=[f'{p}_score' for p in phases])

            cell_cycle_scores["cell_cycle_phase"] = [phases[x] if cell_cycle_scores.iloc[y, x] > 0 else "G1" for y,x in enumerate(np.argmax(cell_cycle_scores.to_numpy(), axis=1))]
            
            self.cell_cycle_scores = cell_cycle_scores
            
            # Save cell cycle scores to file
            out_path = f'{cell_cycle_dir}/cell_cycle_scores.tsv.gz'
            cell_cycle_scores.to_csv(out_path, sep='\t', header=True, index=False)
            
            print('All done!')
            print(f'Cell cycle scores were saved to {out_path}')
            print('########################################\n')

    ### ------------------------------------ ###
    
    def load_cell_cycle_scores(self, data_path):
        
        try:
            
            self.cell_cycle_scores = pd.read_csv(data_path, sep='\t', header=0)
            
        except:
            
            print("ERROR: couldn't read file")
            self.cell_cycle_scores = pd.DataFrame()
    
    ### ------------------------------------ ###

    def score_gene_set(self, gene_set, gene_pool=[], bins=50, ctrl_genes_num=50):

        # Subsetting gene_set and gene_list for detected genes
        gene_set = [gene for gene in gene_set if gene in self.all_genes]
        gene_pool = [gene for gene in gene_pool if gene in self.all_genes]

        # Making sure that gene_set and gene_pool have len > 0
        if not len(gene_set):
            
            print("Error: empty gene_set")
            return None
        
        if not len(gene_pool):
            
            gene_pool = self.all_genes

        # Compute average gene expression levels
        gene_means = np.array(self.all_data.mean(axis=0)).ravel()
        
        # Rank genes based on binned expression level, then for each bin of the genes in the gene_set, pick ctrl_genes random genes for genes with matched binned expression
        bin_size = len(gene_means) // bins
        expression_order = np.argsort(gene_means)
        
        set_order_indexes = []
        for gene in gene_set:
            
            index = expression_order[self.all_genes.index(gene)]
            
            set_order_indexes.append(index)
        
        ctrl_indexes = []
        for index in set_order_indexes:
            
            random_pool = np.where((expression_order >= (index - bin_size // 2)) &
                                   (expression_order <= (index + bin_size // 2)) &
                                   ~(np.isin(expression_order, set_order_indexes)))[0].tolist()
            
            np.random.shuffle(random_pool)
            ctrl_indexes.extend(random_pool[:ctrl_genes_num])
        ctrl_indexes = list(set(ctrl_indexes)) # Removing duplicates

        # Computing the mean of scaled expression values of gene_set genes for each cell
        set_indexes = [self.all_genes.index(g) for g in gene_set if g in self.all_genes]
        set_means = self.scale_features(self.all_data[:, set_indexes]).mean(axis=1)
        #if self.all_data.shape[0] <= self.max_cells_in_memory:
        #    
        #    set_means = self.scale_features(self.all_data[:, set_indexes]).mean(axis=1)
        #
        #else:
        #
        #    set_means = np.array([self.scale_features(self.all_data[:, g]).reshape(-1) for g in set_indexes]).T.mean(axis=1)

        # Computing the mean of scaled expression values of  ctrl_set genes for each cell
        ctrl_means = self.scale_features(self.all_data[:, ctrl_indexes]).mean(axis=1)
        #if self.all_data.shape[0] <= self.max_cells_in_memory:
        #    
        #    ctrl_means = self.scale_features(self.all_data[:, ctrl_indexes]).mean(axis=1)
        #
        #else:
        #
        #    ctrl_means = np.array([self.scale_features(self.all_data[:, g]).reshape(-1) for g in ctrl_indexes]).T.mean(axis=1)
        
        set_score = set_means - ctrl_means
        
        return set_score
    
    ### ------------------------------------ ###
    ### PLOTS                                ###
    ### ------------------------------------ ###
    
    def check_plot_dir(self):
        
        if not hasattr(self, 'plot_dir'):
            
            toggle = 1
        
        elif not exists(self.plot_dir) or not isdir(self.plot_dir):
            
            toggle = 1
        
        else:
            
            toggle = 0
            
        if toggle == 1:
            
            # Create dir to store plots
            plot_dir, count = 'data_plots', 0
            while plot_dir in listdir():
                
                plot_dir = f'data_plots_{count}'
                count += 1
            
            mkdir(plot_dir)
            self.plot_dir = plot_dir
    
    ### ------------------------------------ ###
    
    def plot_cell_cycle(self, datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
        
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='CellCyclePhase', hue_order=['G1PM', 'G1', 'S', 'G2M'], palette='Set1', marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Cell Cycle Phase', markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{self.plot_dir}/Cell_Cycle_Phase_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_cell_type(self, use_subtype=False, cell_types=[], datasets=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]

        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Cell type', hue_order=cell_types, palette=color_palette, marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Cell type', ncol=legend_cols, markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(plot_name, bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_clusters(self, datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
    
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
    
    def plot_datasets(self, datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
    
        # Plotting
        plt.figure(figsize=(5, 5))
        seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Datasets', palette=color_palette, marker='.', s=dot_size, linewidth=0)
        legend = plt.legend(bbox_to_anchor=(1, 1), loc='best', title='Datasets', ncol=legend_cols, markerscale=markerscale)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.savefig(f'{self.plot_dir}/Datasets_UMAP.png', bbox_extra_artists=(legend,), bbox_inches='tight', dpi=300)
        plt.close()
    
    ### ------------------------------------ ###
    
    def plot_gene_expression(self, target_gene, datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
                
            expression = self.all_data[:, gene_index].A.reshape(-1)
        
        if toggle:
            
            # Adding expression data
            plot_data = pd.DataFrame(umap_data, columns=['UMAP_1', 'UMAP_2'])
            plot_data['Expression'] = expression
            
            # Subsetting cells
            plot_data = plot_data.loc[cell_filter,]
            
            # Sort values
            plot_data.sort_values(by="Expression", axis = 0, ascending = True, inplace = True)
            
            # Plotting
            plt.figure(figsize=(6, 5))
            
            if plot_data.Expression.min() != plot_data.Expression.max():
                
                ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Expression', hue_norm=(0, plot_data.Expression.max()), palette='viridis', marker='.', s=dot_size, linewidth=0)
                norm = plt.Normalize(0, plot_data.Expression.max())
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
            
            print(f'Target gene {target_gene}, not found')
    
    ### ------------------------------------ ###
    
    def plot_gene_set_score(self, set_score, gene_set_name='', datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
        
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
    
    def plot_marker_set(self, cell_type_markers, prefix='markers', fig_x_size=10, fig_y_size=6, width_ratios=[6, 0.5, 4]):

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
            
            cl_mean = np.mean(self.all_data[cl_cells,][:, all_markers_idx].toarray(), axis=0).tolist()
            
            plot_data.loc[:, f'cluster_{cl}'] = cl_mean.copy()

        # Remove genes that are not expressed in any sample
        annot_data = annot_data.loc[plot_data.sum(axis=1) > 0,]
        plot_data = plot_data.loc[plot_data.sum(axis=1) > 0,]

        # Clustering
        cm = seaborn.clustermap(plot_data, z_score=0, center=0, vmin=-3, vmax=3, row_colors=annot_data['color'].values, row_cluster=False, col_cluster=True)
        plt.close()
        col_order = cm.dendrogram_col.reordered_ind
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
            
            axes[0].add_patch(plt.Rectangle(xy=(-0.065, i), width=0.05, height=1, color=color, lw=0,
                                            transform=axes[0].get_yaxis_transform(), clip_on=False))
            axes[0].add_patch(plt.Rectangle(xy=(-0.065, i), width=0.05, height=1, edgecolor='black', facecolor='none', lw=1,
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
    
    def plot_trajectories(self, datasets=[], branches=[], pseudotime=False, dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
    
        # Finding cluster centers for plotting branches
        centers = np.array([plot_data.loc[plot_data.Clusters == c, ["UMAP_1", "UMAP_2"]].median(axis = 0) for c in clusters])
    
        # Plotting
        if pseudotime:
            
            # Sort data
            plot_data.sort_values(by="Pseudotime", axis = 0, ascending = True, inplace = True)
            
            # Plot
            plt.figure(figsize=(6, 5))
            ax = seaborn.scatterplot(data=plot_data, x='UMAP_1', y='UMAP_2', hue='Pseudotime', hue_norm=(0, plot_data.Pseudotime.max()), palette='viridis', marker='.', s=dot_size, linewidth=0)
            norm = plt.Normalize(0, plot_data.Pseudotime.max())
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
        
    def plot_variable(self, variable_values, variable_name='Custom variable', categorical=False, sort_values=False, datasets=[], clusters=[], dot_size=1.5, markerscale=2):
        
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
        plot_data = plot_data.loc[cell_filter,]
        
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
            norm = plt.Normalize(plot_data[variable_name].min(), plot_data[variable_name].max())
            sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
            plt.legend().remove()
            plt.xlabel('UMAP 1')
            plt.ylabel('UMAP 2')
            plt.savefig(plot_name, bbox_inches='tight', dpi=300)
            plt.close()
    
### ------------------MAIN------------------ ###

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

from matplotlib import pyplot as plt
from os import listdir, makedirs, mkdir
from os.path import abspath, isdir, exists
from matplotlib import patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from scipy.sparse import csr_matrix, dok_matrix, load_npz, save_npz, vstack
from scipy.stats import hypergeom, mannwhitneyu, zscore
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from statsmodels.stats.multitest import fdrcorrection
