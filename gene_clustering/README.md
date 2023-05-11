# Graph-based gene clustering in PCA space
## Python module for the clustering of genes in PCA space using Leiden algorithm

/// --------------------------------------- ///

### DEPENDENCIES

Python 3.9.7+ &

	h5pyigraph
	matplotlib.pyplot
	numpy
	pandas
	pickle
	random
	seaborn
	sklearn.decomposition.PCA
	sklearn.neighbors.kneighbors_graph
	umap

/// --------------------------------------- ///

### CLASS INIT

cluster_genes class is initialized with the following parameter:

	gene_expression_matrix
					Gene expression matrix in Pandas data-frame format, with
					gene expression values as rows and cells as columns
					Starting columns reporting the gene id or name (can have
					both columns) can be added with names GeneID and GeneSymbol
					Default = None

/// --------------------------------------- ///

### FUNCTIONS

cluster_genes class has the following functions:

/// ---------- ///

load_data(self, file_dir='./')

	Load a previously analyzed dataset, saved as a HDF5 file.

	Attributes:

	file_dir
					Location of the analysis files
					Default = './'

/// ---------- ///

save_data(self, out_dir)

	Saves the analysis to the specified directory.

	Attributes:

	out_dir
					Name of the output directory

/// ---------- ///

reduce_dimensions(self, min_pca_components=20, max_pca_components=50, neighbors=30)

	Reduces the dimensionality of the dataset.
	PCA is performed first, then the number of most informative features is
	found via the elbow method.
	Lastly, UMAP is computed on the selected PCA features.

	Attributes:

	min_pca_components
					Scalar indicating the minimum number of PCA features to use
					Default = 20

	max_pca_components
					Scalar indicating the maximum number of PCA features to use
					Default = 50

	neighbors
					Integer indicating the number of neighbors to be used for
					generating UMAPs
					Default = 30

/// ---------- ///

cluster_genes(self, n_neighbors=10)

	Cluster genes using the Leiden algorithm on a graph of the k
	nearest-neighbors.

	Attributes:

	n_neighbors
					Number of neighbors to use for constructing the graph
					Default = 10

/// ---------- ///

plot_clusters(self, save_dir='./', gene_ids=[], gene_symbols=[], clusters=[], dot_size=1.5)

	Plots the found clusters in UMAP space.

	Attributes:

	save_dir
					Name of the output directory
					Default = []

	gene_ids
					List of gene ids to use (from self.analysis_data)
					Default = []

	gene_symbols
					List of gene symbols to use (from self.analysis_data)
					Default = []

	clusters
					List of clusters to use (from self.analysis_data)
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5
