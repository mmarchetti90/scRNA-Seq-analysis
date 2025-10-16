# 10 Genomics scRNA-Seq analysis
## Python module for the analysis of 10X scRNASeq experiments processed with CellRanger

/// --------------------------------------- ///

### DEPENDENCIES

Python 3.9.16+ &

	gc
	gzip
	igraph
	matplotlib.pyplot
	numpy
	os (listdir, makedirs, mkdir)
	os.path (abspath, isdir, exists)
	pandas
	pickles
	random
	scanorama
	scipy.sparse (csr_matrix, dok_matrix, load_npz, save_npz, vstack)
	scipy.stats (hypergeom, mannwhitneyu)
	seaborn
	sklearn.decomposition import PCA
	sklearn.neighbors.kneighbors_graph
	statsmodels.stats.multitest.fdrcorrection
	umap

/// --------------------------------------- ///

### NOTES

* Designed for RNA Gene Expression data, but can handle ATAC data as well.

/// --------------------------------------- ///

### CLASS INIT

integrated_analysis class is initialized with the following parameter:

	max_cells_in_memory
					Maximum number of cells to load in full matrices
					Default = 100000
	desired_feature_type
					The type of data to extract from MEX files in case of multiomics data (i.e. Gene Expression + ATAC)
					Acceptable values are: 'Gene Expression' and 'Peaks'
					Default = 'Gene Expression'

/// --------------------------------------- ///

### FUNCTIONS

integrated_analysis class has the following functions:

/// ---------- ///

load_manifest(manifest_path)

	Loads the manifest of input files.

	Attributes:

	manifest_path
					Path to the tsv sample manifest

/// ---------- ///

load_all_matrices_in_path(data_path)

	Loads all the sparse matrices from subfolders in the specified path.

	Attributes:

	data_path
					Location of the directory containing the sparse matrices
					subfolders

/// ---------- ///

load_sparse_matrix(main_dir)

	Loads a sparse matrix.
	If data is from a multiome dataset, then the features are filtered for the
	desired type stated at class initialization.

	Attributes:

	main_dir

					Path to directory containing the sparse matrix files

/// ---------- ///

save_sparse_matrix_to_path(main_path, label, barcodes, features, matrix, save_mtx=True)

	Save a sparse matrix to labelled folder.

	Attributes:

	main_path
					Main path where the sparse matrix subfolder will be
					located

	label
					Label for naming the sparse matrix

	barcodes
					List of cells barcodes

	features
					List of gene names

	matrix
					Scipy csr_sparse matrix of counts

	save_mtx
					Set to True to save the counts also in the 10X Genomics
					format
					Default = True

/// ---------- ///

check_matrix_file(path, cells_num, genes_num, max_cells)

	Static method that parses a 10X Genomics mtx file and finds optimal way to
	break it into chunks of self.max_cells_in_memory.

	Attributes:

	path
					Path to mtx file

	cells_num
					Total number of cells

	genes_num
					Total number of genes

	max_cells
					Max cells to load in full matrices

/// ---------- ///

csr_to_mtx(csr, out_name, max_cells=1000)

	Static method to convert a Scipy sparse matrix to 10X Genomics format.

	Attributes:

	csr
					Scipy csr_sparse matrix

	out_name
					Name of output file

	max_cells
					Maximum number of cells to parse at a time
					Default = 1000

/// ---------- ///

full_to_sparse(full)

	Static method to convert a full matrix to Scipy csr_matrix.

	Attributes:

	full
					Full matrix to convert

/// ---------- ///

merge_datasets(datasets_labels, datasets, cells_list, genes_list, desired_genes=[])

	Static method to merge csr_sparse files given a list of cells, and genes.

	Attributes:

	datasets_labels
					List of labels of the datasets to merge

	datasets
					List of csr_matrix datasets

	cells_list
					List of lists of cell names

	genes_list
					List of lists of gene names

	desired_genes
					List of genes to use
					If empty, all genes in common among datasets will be used
					Default = []

/// ---------- ///

read_barcodes_file(path)

	Static method to read list of cell barcodes from file.

	Attributes:

	path
					Path to barcodes file

/// ---------- ///

read_features_file(path)

	Static method to read list of gene names from file.

	Attributes:

	path
					Path to gene names file

/// ---------- ///

read_matrix_chunk(matrix_path, start, stop, matrix_separator):

	Static method to read a chunk of a 10X Genomics formatted mtx file.

	Attributes:

	matrix_path
					Path to mtx file

	start
					First line of the file to be processed

	stop
					Last line of the file to be processed

	matrix_separator
					Field separator used in the mtx file

/// ---------- ///

normalize_data(mtx, skip_first_proportional_fitting=False)

	Normalizes data using PFlog1pPF method.

	Attributes:

	mtx
					csr_sparse file

	skip_first_proportional_fitting
					Set to true if first round of proportional fitting has to
					be skipped (e.g. if counts are in CPM/RPKM/TPM format)
					Default = False

/// ---------- ///

preprocess_samples(min_cell_raw_counts=500, min_detected_genes=500)

	Main function to preprocess samples from dataset manifest file.

	Attributes:

	min_cell_raw_counts
					Minimum number of total raw counts a cell must have to be
					kept
					Default = 500

	min_detected_genes
					Minimum number of genes a cell must have to be kept
					Default = 500

/// ---------- ///

proportional_fitting(mtx, target_size=-1)

	Static method for proportional fitting of cells in a csr_matrix.
	Will scale each cell library size to the median library size.

	Attributes:

	mtx
					csr_sparse matrix

	target_size
					Target size for correction. If set to -1, the median of all
					library sizes will be used
					Default = -1

/// ---------- ///

remove_cells_with_few_genes(brcs, mtx, min_genes)

	Static method to remove cells with too few detected genes.

	Attributes:

	brcs
					List of cell barcodes

	mtx
					csr_sparse matrix

	min_genes
					Min number of detected genes a cell must have

/// ---------- ///

remove_cells_with_low_counts(brcs, mtx, min_cnts)

	Static method to remove cells with too few counts.

	Attributes:

	brcs
					List of cell barcodes

	mtx
					csr_sparse matrix

	min_cnts
					Min total cell counts a cell must have

/// ---------- ///

find_common_highly_variable_genes(data_path, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000, union_mode=False)

	Finds highly variable genes (HVGs) common to all datasets.
	If union mode is True, all HVGs are kept, not just the common ones.
	Based on the plot of exponential of features mean (x axis) vs log of
	features variance/mean (y axis).

	Attributes:

		data_path
					Path to folder containing sparse matrices.

		x_low_cutoff
					Lower cutoff value for x axis
					Default = 0.1

		x_high_cutof
					Higher cutoff value for x axis
					Default = 8

		y_low_cutoff
					Lower cutoff value for y axis
					Default = 1

		y_high_cutoff
					Higher cutoff value for y axis
					Default = 1e9

		nbins
					Number of bins the x axis will be divided into
					Default = 20

		max_features
					Maximum number of variable features saved from each
					dataset
					Default = 3000

		union_mode
					Set to True if all HVGs from all datasets are to be kept,
					not just common ones
					Default = False

/// ---------- ///

find_highly_variable_genes_in_whole_dataset(data_path, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000)

	Finds highly variable genes (HVGs) after merging the datasets.
	Based on the plot of exponential of features mean (x axis) vs log of
	features variance/mean (y axis).

	Attributes:

		data_path
					Path to folder containing sparse matrices.

		x_low_cutoff
					Lower cutoff value for x axis
					Default = 0.1

		x_high_cutof
					Higher cutoff value for x axis
					Default = 8

		y_low_cutoff
					Lower cutoff value for y axis
					Default = 1

		y_high_cutoff
					Higher cutoff value for y axis
					Default = 1e9

		nbins
					Number of bins the x axis will be divided into
					Default = 20

		max_features
					Maximum number of variable features to be saved
					Default = 3000

/// ---------- ///

load_hvgs(data_path)

	Loads a list of HVGs from file.

	Attributes:

	data_path
					Path to HVGs file

/// ---------- ///

find_highly_variable_genes(mtx, x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000)

	Static method to find highly variable genes (HVGs) from a csr_matrix.
	Based on the plot of exponential of features mean (x axis) vs log of
	features variance/mean (y axis).

	Attributes:

		mtx
					csr_sparse matrix

		x_low_cutoff
					Lower cutoff value for x axis
					Default = 0.1

		x_high_cutof
					Higher cutoff value for x axis
					Default = 8

		y_low_cutoff
					Lower cutoff value for y axis
					Default = 1

		y_high_cutoff
					Higher cutoff value for y axis
					Default = 1e9

		nbins
					Number of bins the x axis will be divided into
					Default = 20

		max_features
					Maximum number of variable features to be saved
					Default = 3000

/// ---------- ///

integrate_data(data_path='', use_highly_variable_genes=True, union=False, adjust_library_size=True, use_scanorama_embeddings=True, max_pca_components=50)

	Datasets integration using Scanorama.

	Attributes:

	data_path
					Path to folder containing sparse matrices
					If '', the self.preprocessed_dir will be used
					Default = ''

	use_highly_variable_genes
					If True, HVGs will be used to filter matrices before
					integration (recommended)
					Default = True

	union
					If True, Scanorama will use all genes from all datasets
					(not recommended)
					Default = False

	adjust_library_size
					If True, library sizes will be corrected across datasets
					prior to running Scanorama
					Default = True

	use_scanorama_embeddings
					If True, Scanorama PCA embedding will be used for UMAP
					(recommended)
					Default = True

	max_pca_components
					Default = 50

/// ---------- ///

merge_integrated_datasets(data_path, adjust_library_size=False)

	Merges datasets after integration.
	Called by integrate_data.

	Attributes:

	data_path
					Path to folder containing integrated sparse matrices
					If '', the self.integrated_data_dir will be used
					Default = ''

	adjust_library_size
					If True, library sizes will be corrected across datasets
					Default = False

/// ---------- ///

load_merged_datasets(data_path='')

	Loads a merged integrated dataset.

	Attributes:

	data_path
					Path to folder containing the integrated dataset

/// ---------- ///

reduce_dimensions(use_integrated_pca=True, integrated_pca_path='', use_highly_variable_genes=True, max_pca_components=50, neighbors=30, n_permutations=3)

	Wrapper for umap_embedding and pca_umap_embdedding functions.

	Attributes:

	use_integrated_pca
					If True, the PCA embedding from Scanorama will be used
					Default = True

	integrated_pca_path
	 				Path to Scanorama PCA embdeddings file
	 				Default = ''

	use_highly_variable_genes
	 				If True, use HVGs for PCA embdedding
	 				Default = True

	max_pca_components
	 				Either a scalar indicating the number of PCA features to
					compute or a float in the range 0-1 indicating how much
					variance should be conserved
	 				Default = 50

	neighbors
	 				Integer indicating the number of neighbors to be used for
					generating the UMAP
	 				Default = 30

	n_permutations
	 				If the number of cells in the dataset exceeds
	 				self.max_cells_in_memory, PCA and UMAP will be fitted using
	 				batches of cells, with all cells used n_permutations times
	 				Default = 3

/// ---------- ///

pca_umap_embdedding(use_highly_variable_genes, max_pca_components, neighbors, n_permutations)

	Reduces the dimensionality of the dataset.
	PCA is performed first, then the number of most informative features is
	found via the elbow method.
	Lastly, UMAP is computed on the selected PCA features.

	Attributes:

	use_highly_variable_genes
					If True, only the most variable features are used for PCA

	max_pca_components
					Either a scalar indicating the number of PCA features to
					compute or a float in the range 0-1 indicating how much
					variance should be conserved

	neighbors
					Integer indicating the number of neighbors to be used for
					generating the UMAP

	n_permutations
	 				If the number of cells in the dataset exceeds
	 				self.max_cells_in_memory, PCA and UMAP will be fitted using
	 				batches of cells, with all cells used n_permutations times

/// ---------- ///

scale_features(mtx, gene_mean=[], gene_std=[])

	Scales features by mean centering first then dividing by standard deviation
	Required for PCA step.

	Attributes:

	mtx
					csr_sparse matrix

	gene_mean
					Pre-calculated gene means
					If [], new values will be calculated using available cells
					Default = []

	gene_std
					Pre-calculated gene standard deviations
					If [], new values will be calculated using available cells
					Default = []

/// ---------- ///

umap_embedding(self, integrated_pca_path, neighbors, n_permutations)

	Reduces the dimensionality of the dataset using UMAP on Scanorama PCA
	embeddings.

	Attributes:

	integrated_pca_path
	 				Path to Scanorama PCA embdeddings file
	 				Default = ''

	neighbors
	 				Integer indicating the number of neighbors to be used for
					generating the UMAP

	n_permutations
	 				If the number of cells in the dataset exceeds
	 				self.max_cells_in_memory, PCA and UMAP will be fitted using
	 				batches of cells, with all cells used n_permutations times

/// ---------- ///

get_gene_mean_and_std(mtx)

	Static method to calculate gene means and standard deviations from a
	csr_sparse matrix.

	Attributes:

	mtx
					csr_sparse matrix

/// ---------- ///

load_reduced_dimensions(lower_dimensions_dir)

	Static method to load PCA and UMAP embdeddings.

	Attributes:

	lower_dimensions_dir
					Path to directory contaning PCA and UMAP data

/// ---------- ///

make_cell_subsets(N, max_cells)

	Static method to create random groups of max_cells.

	Attributes:

	N
					Total number of cells available

	max_cells
					Maximum number of cells in a group

/// ---------- ///

check_cluster_composition()

	Checks the distribution of cells within clusters using a hypergeometric
	test.
	If any dataset is significantly enriched in a cluster, a flag is raised.

	Attributes:

	

/// ---------- ///

cluster_cells(self, n_neighbors=10, resolution_parameter=1.0, beta=0.01)

	Cluster cells using the Leiden algorithm on a graph of the k
	nearest-neighbors.

	Attributes:

	n_neighbors
					Number of neighbors to use for constructing the graph
					Default = 10

	resolution_parameter
					Higher resolutions lead to more smaller communities, while lower resolutions lead to fewer larger communities
					Default = 1.0

	beta
					Parameter affecting the randomness in the Leiden algorithm
					Default = 0.01

/// ---------- ///

find_clusters_markers(self, min_fc=0.25, min_pct=0.1, min_cells_in_cluster=50)

	Finds the features differentially expressed in each cluster.

	Attributes:

	min_fc
					Minimum Log2 fold-change between the two cell populations
					for the feature to be considered
					Default = 0.25

	min_pct
					Minimum percentage of cells in either population expressing
					the feature for it to be considered
					Default = 0.1

/// ---------- ///

load_clusters(self, data_path)

	Loads clustering information from file.

	Attributes:

	data_path
					Path to file containing cluster embeddings

/// ---------- ///

load_cluster_identities(self, data_path)

	Loads cluster identities from file.

	Attributes:

	data_path
					Path to file containing cluster identities

/// ---------- ///

load_cluster_markers(self, data_path)

	Loads cluster markers from file.

	Attributes:

	data_path
					Path to file containing cluster markers

/// ---------- ///

score_cluster_identity(self, ref_expression, ref_samples, clusters_id=[])

	Scores the identity of clusters by comparing its markers (found with
	find_clusters_markers) to a reference expression dataset of annotated
	cells.

	Attributes

	ref_expression
					Reference expression data frame
					Must have the following columns:

						"gene_name", i.e. gene symbols or identifiers matching
						the self.cluster_markers.gene_name

						A column for each cell type with, for each feature, the
						log2 of the ratio between feature expression in that
						cell type and the mean expression across all cell types

	ref_samples
					Reference samples data frame
					Must have the following columns:

						"reference_cell_type", i.e. main cell type name

						"long_name", i.e. the name of the specific cell
						(sub)type

	clusters_id
					List of clusters to use
					Default = []

/// ---------- ///

compare_populations(gene_names, data, pop1, pop2, max_cells_in_memory=10000, min_fc=0.25, min_pct=0.1)

	Compares the expression of features between cells in one cluster and cells
	from every other cluster using a Mann-Whitney U Test with BH correction.

	Attributes:

	pop1, pop2
					List of cell names used to subset data

	min_fc
					Minimum Log2 fold-change between the two cell populations
					for the feature to be considered
					Default = 0.25

	min_pct
					Minimum percentage of cells in either population expressing
					the feature for it to be considered
					Default = 0.1

/// ---------- ///

create_mst(self, clusters_id=[])

	Trajectory analysis using clusters centroids in PCA space to create a
	minimum spanning tree with Kruskal algorithm.

	Attributes:

	clusters
					List of clusters to use
					Default = []

/// ---------- ///

generate_pseudotime(self)

	Pseudotime calculation based on the distance from the root vertex of the
	MST tree.

/// ---------- ///

load_trajectory_branches(self, data_path)

	Loads MST branched generated by create_mst.

	Attributes:

	data_path
					Path to file containing MST branches

/// ---------- ///

load_trajectory_pseudotime(self, data_path)

	Loads cells pseudotime values calculated by generate_pseudotime.

	Attributes:

	data_path
					Path to file containing cells pseudotime values

/// ---------- ///

score_cell_cycle(self, g1pm_features, s_features, g2m_features, gene_pool=[], bins=50, ctrl_genes_num=50)

	Wrapper for score_gene_set that scores G1 postmitotic, G1S, and G2M
	phase-specific genes, similarly to Seurat.
	Gene set with the highest positive value is used to assign a predicted
	cell cycle phase. If all gene sets have a negative score, cell is
	assigned G1 phase.

	Attributes:

	g1pm_features, s_features, g2m_features
					Sets of gene IDs or symbols specific for G1 postmitotic
					cells, S-phase cells, and G2M-phase cells

	gene_pool
					Set of genes to use as random control
					If not specified, the whole feature space is used
					Default = []

	bins
					Number specifying into how many sets the range of average
					gene expression values should be divided into
					For each target feature belonging to the gene set analyzed,
					ctrl_genes_num genes will be randomly selected as controls
					of matched expression
					Default = 50

	ctrl_genes_num
					Number of control genes to pick from each expression bin
					Default = 50

/// ---------- ///

load_cell_cycle_scores(self, data_path)

	Loads cell cycle scores generated by score_cell_cycle.

	Attributes:

	data_path
					Path to cell cycle scores file

/// ---------- ///

score_gene_set(self, gene_set, gene_pool=[], bins=50, ctrl_genes_num=50)

	Scoring expression of a set of genes versus a random set of genes with
	similar expression range in the whole cell population, similarly to Seurat.

	Attributes:

	gene_set
					Set of gene IDs or symbols of interest

	gene_pool
					Set of genes to use as random control
					If not specified, the whole feature space is used
					Default = []

	bins
					Number specifying into how many sets the range of average
					gene expression values should be divided into
					For each target feature belonging to the gene set analyzed,
					ctrl_genes_num genes will be randomly selected as controls
					of matched expression
					Default = 50

	ctrl_genes_num
					Number of control genes to pick from each expression bin
					Default = 50

/// ---------- ///

check_plot_dir(self)

	Checks if the parameter self.plot_dir exists and creates a directory for
	storing plots, otherwise.

/// ---------- ///

plot_cell_cycle(self, datasets=[], clusters=[], dot_size=1.5)

	Plots the predicted cell cycle phase generated by score_cell_cycle in UMAP
	space.

	Attributes:

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_cell_type(self, use_subtype=False, datasets=[], cell_types=[], dot_size=1.5)

	Plots the inferred cluster identity generated by score_cluster_identity in
	UMAP space.

	Attributes:

	use_subtype
					Set to true to use cell sub-types instead of main types.
					Default = False

	datasets
					List of labels of datasets to plot
					Default = []

	cell_types
					List of cell (sub)types to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_clusters(self, datasets=[], clusters=[], dot_size=1.5)

	Plots the found clusters in UMAP space.

	Attributes:

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_datasets(self, datasets=[], clusters=[], dot_size=1.5)

	Plots the datasets distribution in UMAP space.

	Attributes:

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_gene_expression(self, target_gene, datasets=[], clusters=[], dot_size=1.5)

	Plots the expression of a desired target gene in UMAP space.

	Attributes:

	target_gene
					Feature identifier

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_gene_set_score(self, set_score, gene_set_name='', datasets=[], clusters=[], dot_size=1.5)

	Plots the expression of a gene set in UMAP space.

	Attributes:

	set_score
					Expression score of the gene set in each cell, generated by
					score_gene_set.

	gene_set_name
					Name of the gene set
					Default = ''

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	dot_size
					Size of the scatterplor dots
					Default = 1.5

/// ---------- ///

plot_marker_set(self, cell_type_markers, prefix='markers', fig_x_size=10, fig_y_size=6, width_ratios=[6, 0.5, 4])

	Plots a heatmap of a desired set of markers, annotated for the corresponding cell type.

	Attributes:

	cell_type_markers
					Pandas DataFrame with two columns: cell_name, and markers.

	prefix
					Plot prefix.
					Default = 'markers'

	fig_x_size
					Figure width
					Default = 10

	fig_y_size
					Figure height
					Default = 6

	width_ratio
					Width ratios of subplots
					Default = [6, 0.5, 4]

/// ---------- ///

plot_trajectories(self, datasets=[], clusters=[], pseudotime=False, dot_size=1.5)

	Plots the trajectories from create_mst in UMAP space.
	Only the clusters included in the MST will be plotted.

	Attributes:

	set_score
					Expression score of the gene set in each cell, generated by
					score_gene_set.

	datasets
					List of labels of datasets to plot
					Default = []

	clusters
					List of clusters to plot
					Default = []

	pseudotime
					If True, will use pseudotime for coloring individual cells
					Default = False

	dot_size
					Size of the scatterplor dots
					Default = 1.5
