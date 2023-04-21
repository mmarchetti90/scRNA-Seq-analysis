# 10 Genomics scRNA-Seq analysis
## Python module for the analysis of 10X scRNASeq experiments processed with CellRanger

/// --------------------------------------- ///

### DEPENDENCIES

Python 3.9.7+ &

	h5pyigraph
	matplotlib.pyplot
	numpy
	os.listdir
	pandas
	pickle
	random
	seaborn
	scipy.stats.mannwhitneyu
	sklearn.decomposition.PCA
	sklearn.neighbors.kneighbors_graph
	statsmodels.stats.multitest.fdrcorrection
	umap

/// --------------------------------------- ///

### FUNCTIONS

tenx_scranseq class has the following functions:

/// ---------- ///

load_processed_data(file_name, file_dir='./')

	Load a previously analyzed dataset, saved as a HDF5 file.

	Attributes:

	file_name
					Name of the HDF5 file

	file_dir
					Location of the HDF5 file, as well as Pickle files where
					the PCA and UMAP models are stored.
					Default = './'

/// ---------- ///

saveClassObject(object_name="scRNASeq_Analysis")

	Saves the analysis.
	Outputs a HDF5 file containing all the info from the class instance, as
	well as Pickle files for the PCA and UMAP models.

	Attributes:

	object_name
					Name of the output HDF5 file
					Default = "scRNASeq_Analysis"

/// ---------- ///

import_raw_data(input_dir="./", min_counts=500)

	Import raw counts from CellRanger.

	Attributes:

		input_dir
					Location of the HDF5 files
					Default = "./"

		min_counts
					Minimum total count for individual cells
					Default = 500

/// ---------- ///

recreate_full_matrix(file_id, file, min_counts)

	Static method to read raw counts HDF5 files from CellRanger and returns a
	full matrix.

	Attributes:

		file_id	
					Name of the HDF5 file to be used for cells' metadata

		file
					Path to the HDF5 file

		min_counts
					Minimum total count for individual cells

/// ---------- ///

normalize_counts(scale_factor=10000)

	LogNormalizes raw counts, similar to Seurat's NormalizeData.

	Attributes:

		scale_factor
					Factor for library size scaling
					Default = 10000

/// ---------- ///

find_most_variable_features(x_low_cutoff=0.1, x_high_cutof=8, y_low_cutoff=1, y_high_cutoff=1e9, nbins=20, max_features=3000)

	Finds most variable genes in the dataset, similarly to Seurat.
	Based on the plot of exponential of features mean (x axis) vs log of
	features variance/mean (y axis).

	Attributes:

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
					Maximum number of variable features saved
					Default = 3000

/// ---------- ///

scale_features()

	Scales features by mean centering first then dividing by standard deviation.
	Required for PCA step.

/// ---------- ///

reduce_dimensions(use_variable_features=True, pca_components=50, neighbors=[5, 30, 50])

	Reduces the dimensionality of the dataset.
	PCA is performed first, then the number of most informative features is
	found via the elbow method.
	Lastly, UMAP is computed on the selected PCA features.

	Attributes:

	use_variable_features
					If True, only the most variable features are used for PCA
					Default = True

	pca_components
					Either a scalar indicating the number of PCA features to
					compute or a float in the range 0-1 indicating how much
					variance should be conserved
					Default = 50

	neighbors
					List of integers indicating the number of neighbors to be
					used for generating UMAPs
					A UMAP will be generated for each value in the list
					Default = [5, 30, 50]

/// ---------- ///

cluster_cells(n_neighbors=10)

	Cluster cells using the Leiden algorithm on a graph of the k
	nearest-neighbors.

	Attributes:

	n_neighbors
					Number of neighbors to use for constructing the graph
					Default = 10

/// ---------- ///

find_clusters_markers(min_fc=0.25, min_pct=0.1)

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

compare_populations(pop1, pop2, min_fc=0.25, min_pct=0.1)

	Compares the expression of features between cells in one cluster and cells
	from every other cluster using a Mann-Whitney U Test with BH correction.

	Attributes:

	pop1, pop2
					List of cell names used to subset self.normalized_counts

	min_fc
					Minimum Log2 fold-change between the two cell populations
					for the feature to be considered
					Default = 0.25

	min_pct
					Minimum percentage of cells in either population expressing
					the feature for it to be considered
					Default = 0.1

/// ---------- ///

score_cluster_identity(ref_expression, ref_samples, clusters=[])

	Scores the identity of clusters by comparing its markers (found with
	find_clusters_markers) to a reference expression dataset of annotated
	cells.

	Attributes:

	ref_expression
					Reference expression data frame
					Must have the following columns:

						"GeneName", i.e. gene symbols or identifiers matching
						the self.cluster_markers.GeneName

						A column for each cell type with, for each feature, the
						log2 of the ratio between feature expression in that
						cell type and the mean expression across all cell types

	ref_samples
					Reference samples data frame
					Must have the following columns:

						"reference_cell_type", i.e. main cell type name

						"short_name", i.e. the abbreviated name of the specific
						cell (sub)type

						"long_name", i.e. the name of the specific cell
						(sub)type

	clusters
					List of clusters to use
					Default = []

/// ---------- ///

create_mst(clusters=[])

	Trajectory analysis using clusters centroids in PCA space to create a
	minimum spanning tree with Kruskal algorithm.

	Attributes:

	clusters
					List of clusters to use
					Default = []

/// ---------- ///

generate_pseudotime()

	Pseudotime calculation based on the distance from the root vertex of the
	MST tree.

/// ---------- ///

score_cell_cycle(g1pm_features, s_features, g2m_features, gene_pool=[], bins=50, ctrl_genes_num=50)

	Wrapper for score_gene_set that scores G1 postmitotic, G1S, and G2M
	phase-specific genes, similarly to Seurat.
	Gene set with the highest positive value is used to assign a predicted
	cell cycle phase. If all gene sets have a negative score, cell is
	assigned G1 phase.

	Attributes:

	g1pm_features, s_features, g2m_features
					Sets of gene IDs or symbols specific for G1 postmitotic
					cells, S-phase cells, and G2M-phase cells.

	gene_pool
					Set of genes to use as random control.
					If not specified, the whole feature space is used.
					Default = []

	bins
					Number specifying into how many sets the range of average
					gene expression values should be divided into.
					For each target feature belonging to the gene set analyzed,
					ctrl_genes_num genes will be randomly selected as controls
					of matched expression.
					Default = 50

	ctrl_genes_num
					Number of control genes to pick from each expression bin.
					Default = 50

/// ---------- ///

score_gene_set(gene_set, gene_pool=[], bins=50, ctrl_genes_num=50)

	Scoring expression of a set of genes versus a random set of genes with
	similar expression range in the whole cell population, similarly to Seurat.

	Attributes:

	gene_set
					Set of gene IDs or symbols of interest.

	gene_pool
					Set of genes to use as random control.
					If not specified, the whole feature space is used.
					Default = []

	bins
					Number specifying into how many sets the range of average
					gene expression values should be divided into.
					For each target feature belonging to the gene set analyzed,
					ctrl_genes_num genes will be randomly selected as controls
					of matched expression.
					Default = 50
			   
	ctrl_genes_num
					Number of control genes to pick from each expression bin.
					Default = 50

/// ---------- ///

plot_gene_expression(target_gene, cells=[], samples=[], groups=[], clusters=[], dot_size=1.5)

	Plots the expression of a desired target gene in UMAP space.

	Attributes:

	target_gene
					Feature identifier (from either self.raw_counts.GeneID or
					self.raw_counts.GeneID)

	cells
					List of cell names used to use (from self.metadata)
					Default = []

	samples
					List of samples identifiers to use (from self.metadata)
					Default = []

	groups
					List of samples groups identifiers to use
					(from self.metadata)
					Default = []

	clusters
					List of clusters to use (from self.metadata)
					Default = []

	dot_size
					Size of the scatterplor dots.
					Default = 1.5

/// ---------- ///

plot_clusters(cells=[], samples=[], groups=[], clusters=[], dot_size=1.5)

	Plots the found clusters in UMAP space.

	Attributes:

	cells
					List of cell names used to use (from self.metadata)
					Default = []

	samples
					List of samples identifiers to use (from self.metadata)
					Default = []

	groups
					List of samples groups identifiers to use
					(from self.metadata)
					Default = []

	clusters
					List of clusters to use (from self.metadata)
					Default = []

	dot_size
					Size of the scatterplor dots.
					Default = 1.5

/// ---------- ///

plot_trajectories(cells=[], samples=[], groups=[], pseudotime=False, dot_size=1.5)

	Plots the trajectories from create_mst in UMAP space.
	Only the clusters included in the MST will be plotted.

	Attributes:

	cells
					List of cell names used to use (from self.metadata)
					Default = []

	samples
					List of samples identifiers to use (from self.metadata)
					Default = []

	groups
					List of samples groups identifiers to use
					(from self.metadata)
					Default = []

	pseudotime
					If True, will use pseudotime for coloring individual cells
					Default = False

	dot_size
					Size of the scatterplor dots.
					Default = 1.5

/// ---------- ///

plot_gene_set_score(set_score, cells=[], samples=[], groups=[], clusters=[], dot_size=1.5)

	Plots the expression of a gene set in UMAP space.

	Attributes:

	set_score
					Expression score of the gene set in each cell, generated by
					score_gene_set.

	cells
					List of cell names used to use (from self.metadata)
					Default = []

	samples
					List of samples identifiers to use (from self.metadata)
					Default = []

	groups
					List of samples groups identifiers to use
					(from self.metadata)
					Default = []

	dot_size
					Size of the scatterplor dots.
					Default = 1.5

/// ---------- ///

plot_cell_cycle(cells=[], samples=[], groups=[], clusters=[], dot_size=1.5)

	Plots the predicted cell cycle phase generated by score_cell_cycle in UMAP
	space.

	Attributes:

	cells
					List of cell names used to use (from self.metadata)
					Default = []

	samples
					List of samples identifiers to use (from self.metadata)
					Default = []

	groups
					List of samples groups identifiers to use
					(from self.metadata)
					Default = []

	dot_size
					Size of the scatterplor dots.
					Default = 1.5
