# Single-cell RNASeq of diverse datasets

## 0 - Intro

This tutorial will guide you in how to integrate and analyze datasets from different origins.
The pipeline is written entirely in Python and requires several packages (see README file for details), but can also be run using a docker image (see dockerfile). Unfortunately, the Umap package in the docker image is not working properly, so any steps requiring Umap dimensionality reduction (i.e. mostly visualizations) cannot be run using the image. However, the data integration, which is the most intensive, does not require the Umap module. So, the docker image can be used on a HPC cluster to avoid installing the required Python packages.

### 0.1 - Tutorial datasets

This tutorial uses two published and publicly available datasets:

1. GSE125188 (Zhao et al 2020), CD45+ cells isolated from liver from 3 patients. Data is raw counts stored in the MEX format used by 10X CellRanger. This is comprised of 3 tables for storing barcodes, features (e.g. genes), and counts as sparse matrices.

2. GSE115469 (MacParland et al 2018), total liver homogenate from 5 patients. Data is pooling normalized counts and stored in a large csv table.

Both datasets should be downloaded from GEO and stored locally as follows:
For GSE125188, download the barcodes, genes, and matrix files and store them in a newly created GSE125188_cd45+ directory.
For GSE115469, use the dataset_reformatting_scripts/reformat_GSE115469.py script to convert it to a format usable by the pipeline.

### 0.2 - Datasets formatting

The pipeline can processed data stored in two ways:

1. Tab-separated matrices of counts, where the first column must be named 'GeneSymbol' and contain (you guessed it) gene symbols. The following columns contain counts for individual cells. Can be compressed with gzip.

2. MEX format used by 10X CellRanger, with 3 needed files for barcodes, features (e.g. genes), and counts, which must contain the keywords 'barcodes', 'features', and 'matrix' in their names, respectively. Can be compressed with gzip.

### 0.3 - Pipeline input

The pipeline requires a datasets manifest as input, which specifies several fields of which 3 are required:

- data_label = label with which cells from a specific dataset will be named.

- dataset_path = full path to the dataset (see below for details).

- normalization = type of normalization used on the data. Can be any of:
    - "none" = raw counts
    - "cpm" = counts per million
    - "rpkm" = reads per kilobase per million mapped reads
    - "pool_norm" = pooling normalization
    - "pflog1ppf" = proportional fitting followed by log1p followed by proportional fitting

Additional fields can be added, but will be ignored by the pipeline. See attached example for details.

Paths for tab-separated matrices are full paths to the matrix itself, while paths for MEX formatted data must direct to the folder containing the barcodes, features, and counts tables.


## 1 - Data preprocessing

### 1.1 - Overview

Typical data processing always includes the following several steps:

- Data normalization and filtering

- Highly variable genes detection

- Datasets integration

- Dimensionality reduction

- Cell clustering

### 1.2 - Setup

First of all, import the package:

```
import integrated_scrnaseq_analysis as isa
```

Next, we need to initialize an integrated_analysis class, specifying how many cells we want to keep in memory at most. If your machine is struggling with large datasets, lowering this number can help reduce the memory footprint. Keep in mind that the pipeline stores count data in sparse matrices from the scipy package, which are quite memory efficient, so even a modest laptop can handle 100000 cells at a time.
After initializing the class, the manifest is loaded.

```
analysis = isa.integrated_analysis(max_cells_in_memory=100000)

analysis.load_manifest(manifest_path)
```

### 1.3 - Sample preprocessing

In the first step of data preprocessing, for each dataset, cells with too few detected genes are discarded. Then data is normalized, except when a sample was already normalized using 'pflog1ppf' or 'pool_norm' methods. Note that if raw count data is available, then the min_cell_raw_counts parameter is also used to discard cells with too few reads. Lastly, preprocessed datasets are saved locally as sparse matrices in a newly created directory with prefix '1_preprocessed_data'. The path to this directory is stored as a string in the analysis.preprocessed_dir parameter.
All these steps are wrapped in the following function.

```
analysis.preprocess_samples(min_cell_raw_counts=500, min_detected_genes=500)
```

Note that this step is quite time-consuming. However, since preprocessed, ready for integration, datasets are saved locally, there's no need to reprocess every dataset when new ones need to be added to the analysis. One could simply preprocess the new datasets, then store them in the same directory as others previously preprocessed, specifying this path in the analysis.preprocessed_dir variable.

### 1.4 - Highly variable genes detection

Prior to any downstream analysis, datasets will have to be integrated to reduce variability across samples (e.g. due to different batches or sample preps). This can be done using all the genes from the datasets, but focusing on genes that are highly variable (HVGs) yields better results.
The find_common_highly_variable_genes function is used to find genes that are highly variable in each individual dataset. If the attribute 'union' is set to True, then all HVGs from individual datasets will be combined, while with the default False value, only common ones across datasets will be used. It is recommended to set 'union' to False to focus on common elements between the datasets, but one should use the union mode if too few genes are found. Alternatively, the find_highly_variable_genes_in_whole_dataset function can be used to find highly variable genes after merging all datasets into one if the find_common_highly_variable_genes returns few HVGs.
The function will generate a list of HVGs stored as a list in the analysis.hvgs and locally in a newly created directory with prefix '2_hvgs'.
The code snippet below also shows how to manually add genes to the analysis.hvgs list, in case some gene of interest is missing from it.

```
analysis.find_common_highly_variable_genes(data_path=analysis.preprocessed_dir, union_mode=True)
analysis.hvgs += ['CCR5', 'CCR6', 'CD3D', 'CD3E', 'CD3G', 'CD4']
```

### 1.5 - Datasets integration

The Python package Scanorama is used for datasets integration. The integrate_data function is a wrapper for Scanorama and creates both PCA embeddings and corrected counts, which will be stored locally in a newly created directory with prefix '3_integrated_data'. Then, data is merged to a single dataset and saved locally in a newly created '4_merged_data' directory.
The data_path attribute is required and specifies the location of the preprocessing directory. Note that this info is also stored in the analysis.preprocessed_dir, so an empty string can also be passed as data_path, in which case the function will use the analysis.preprocessed_dir path.
The 'use_highly_variable_genes' allows the datasets to be filtered for the HVGs detected above prior to integration. Then, to make sure that Scanorama only uses genes detected in all datasets (recommended), we need to set the 'union' attribute to True.
Scanorama will not only correct count matrices, but also generate PCA embedding. It is recommended to use these for downstream analyses (e.g. clustering, UMAP, ...), but by setting the use_scanorama_embeddings attribute to False, one can instead run a PCA at a later time.

```
analysis.integrate_data(data_path='', use_highly_variable_genes=True, union=False, use_scanorama_embeddings=True, max_pca_components=50)
```

### 1.6 - Dimensionality reduction

Dimensionality reduction is another essential step in scRNASeq analysis and is run using the reduce_dimensions wrapper. This function comprises two nested steps: PCA and UMAP. Note that only the former is used by downstream analyses (e.g. clustering and trajectory analysis), while Umaps are only used for visualizations. Also, since dimensionality reduction relies on the Umap package, the Docker image will not be useful for this and next steps requiring Umap.
If the use_integrated_pca attribute is set to True (default), then the PCA embeddings from Scanorama will be used. Otherwise, PCA will be run before Umap, in which case focusing on HVGs is recommended (set the use_highly_variable_genes attribute to True). If using Scanorama PCA embeddings, these are normally taken from the analysis.pca_data variable (derived from the analysis.integrate_data function), but, if this variable is missing (e.g. analysis was interruped and is now being continued), then the path to the file containing the Scanorama PCA embeddings can be specified using the integrated_pca_path attribute.
In the code below, the reduce_dimensions function is called using Scanorama-generated PCA embeddings.

```
analysis.reduce_dimensions(use_integrated_pca=True)

analysis.plot_datasets(dot_size=10)
```

### 1.7 - Cell clustering

A graph-based approach with the Leiden algorithm is used for cell clustering. Note that modifying the n_neighbors attribute will affect the 'granularity' of the clustering. Clusters will be stored in analysis.clusters and locally in a newly created directory with prefix '6_clustering'.
The script below also shows how to find markers for each cluster, i.e. genes differentially expressed in the cells of a specific cluster vs all the other cells. Cluster markers will be stored in analysis.cluster_markers and locally in a newly created directory with prefix '7_cluster_markers'.

```
analysis.cluster_cells(n_neighbors=20)

analysis.plot_clusters()

analysis.find_clusters_markers(min_fc=0.25, min_pct=0.1, min_cells_in_cluster=50)
```


## 2 - Other analyses

### 2.1 - Trajectory analysis

The following code snippet is used to run a Kruskal Minimum Spanning Tree algorithm, creating a minimum spanning tree (i.e. a trajectory) to connect clusters based on their PCA transformed gene expression matrix. In short, a branching tree is generated which connects each cluster while minimizing total branch length. This is useful to have an idea of which cell populations are closer to one another in terms of gene expression. The create_mst function below will create a tree connecting all clusters, but the tree can also be created on a subset of clusters which can be passed as a list to the create_mst as the clusters_id attribute.
The following functions, generate_pseudotime, can be used to visualize the distance of individual cells from the cluster at the base of the tree, giving a more visual idea of 'distance' between cells' expression.
Data generated by these functions is stored in the analysis.branches and analysis.pseudotime variables and locally in a newly created directory with prefix '9_trajectories'.

```
analysis.create_mst()

analysis.generate_pseudotime()

analysis.plot_trajectories(pseudotime=False)
analysis.plot_trajectories(pseudotime=True)
```

### 2.2 - Differential expression between cell populations

The function below is used to compare two cell populations, finding differentially expressed genes. The two populations are defined by passing two filters, i.e. lists of booleans. This is useful, for example, to compare two different populations (e.g. healthy and diseased) within the same cluster.
The example below compared cells from cluster 1 and 2. Indeed, the find_clusters_markers function described above is just a wrapper for compare_populations.
Note that gene names and expressions have to be passed to the function. This information can be retrieved from the analysis.all_genes variables, respectively.

```
pop1 = (analysis.clusters == 1)
pop2 = (analysis.clusters == 2)

analysis.compare_populations(analysis.all_genes, analysis.all_data, pop1, pop2)
```

### 2.3 - Marker signature analysis

When interested in a specific set of genes, the most intuitive approach is to plot each of them individually. The script below shows just that using the plot_gene_expression function.
Another approach, easier to visualize, is to compute a score of the whole gene set. This is expressed for each individual cell as the ratio of the expression of the genes in the set and the expression of a set of random genes that, in the whole cell population, have similar levels to the genes in the set. So, the score indicates wheter the set is highly or lowly expressed in an individual cell compared to the general population. The score_gene_set returns is used for this purpose.

```
mait_gene_set = ['CD3D', 'CD3E', 'CD3G', 'KLRB1', 'TRAV1-2']

for gene in mait_gene_set:

    analysis.plot_gene_expression(gene)

mait_gene_set_scores = analysis.score_gene_set(mait_gene_set)

analysis.plot_gene_set_score(mait_gene_set_scores, gene_set_name='mait_gene_set')
```

### 2.4 - Cell Cycle Analysis

The score_cell_cycle function is a wrapper for the score_gene_set, calculating a score for G1 post-mitotic (G1PM), S, and G2M cell cycle genes (lists of genes have to be provided for each category). Then, the highest score is used to define the cell cycle phase for each cell.

```
analysis.score_cell_cycle(g1pm_genes, s_genes, g2m_genes)

analysis.plot_cell_cycle()
```

### 2.5 - Cell type identification

When searching for rare cell types, different approaches may yield different, possibly contradicting results. The find_specific_cell_type module can be used to combine different metrics to get a more robust prediction. The find_cell_type class can use both methods based on the signature of a selected gene set (AUC and gene set score) and methods based on a reference (euclidean distance and Pearson's and Spearman's correlation). Metrics are computed for each cluster, normalized to a 0-1 scale via softmax (or softmin for euclidean distance), then compounded for each cluster individually.

Continuing with the integrated_analysis class analyzed above ('analysis'), we first need to load the find_specific_cell_type module, the gene set, and the reference expression data.

```
import find_specific_cell_type as fsct
import pandas as pd

# Mait min signature
mait_gene_set = ['CD3D', 'CD3E', 'CD3G', 'KLRB1', 'TRAV1-2']

# Loading MAIT reference data
ref_expr = pd.read_csv('mait_sc_reference.tsv.gz', sep='\t')
ref_expr = ref_expr.loc[ref_expr.GeneName.isin(analysis.all_genes), ] # Filter ref_expr, keeping only genes in the 'analysis' object
```

Then, the find_cell_type object is initialized using the counts, gene names, and clusters from the 'analysis' object.

```
mait_search = fsct.find_cell_type(analysis.all_data, analysis.all_genes, analysis.clusters, mait_gene_set, ref_expr)
```

The analyze_data function is a wrapper to run the individual analysis.

```
results = mait_search.analyze_data(max_cells_in_memory=200000)
results.to_csv('mait_search.tsv', sep='\t', index=False)
```

Lastly, the diagnostic_plots function can be invoked to plot useful visualizations of the individual and compounded metrics.

```
mait_search.diagnostic_plots()
```

## 3 - Data subsetting

The data formats used by the pipeline make it easy to subset a population of cells to create a new dataset. This section will show you how to extract cells from cluster 1 of the 'analysis' object described in the sections above and create a new dataset from it.

First of all, we need to load a few more modules.

```
import gzip
import numpy as np
import pandas as pd
from scipy.sparse import save_npz
```

Now we initialize a new integrated_analysis class.

```
analysis_subset = isa.integrated_analysis(max_cells_in_memory=100000)
```

Then we create a boolean filter to select cells from cluster 1.

```
cell_filter = (analysis.clusters == 1)
```

Now we can subset cells and their data using this filter.

```
all_data_subset = analysis.all_data[cell_filter,].copy()
all_cells_subset = list(np.array(analysis.all_cells)[cell_filter])
```

We'll also have to create a new datasets_metadata object, which stores info about the datasets that make up the integrated data.

```
datasets_metadata_subset = {'datasets_labels' : [], 'cell_num' : []}
for dataset_label in analysis.datasets_metadata.datasets_labels:
    
    cell_num = sum([1 for cell in all_cells_subset if cell.endswith(dataset_label)])
    
    if cell_num > 0:
        
        datasets_metadata_subset['datasets_labels'].append(dataset_label)
        datasets_metadata_subset['cell_num'].append(cell_num)

datasets_metadata_subset = isa.pd.DataFrame(datasets_metadata_subset)
```

The analysis_subset object can now be filled with the subset data above.

```
analysis_subset.all_data = all_data_subset.copy()
analysis_subset.all_cells = all_cells_subset.copy()
analysis_subset.all_genes = analysis.all_genes.copy() # Genes were not subset, so will just be copied
analysis_subset.datasets_metadata = datasets_metadata_subset.copy()
```

We'll also need to create local copies so that we'll be able to reload the object at a later time (see next below).

```
# Save metadata
datasets_metadata_subset.to_csv('datasets_metadata.tsv.gz', sep='\t', index=None)

# Save barcodes
out_text = ('\n'.join(all_cells_subset)).encode()
with gzip.open('all_barcodes.tsv.gz', 'wb') as output:
    
    output.write(out_text)

# Save features
out_text = ('\n'.join(analysis.all_genes)).encode()
with gzip.open('all_features.tsv.gz', 'wb') as output:
    
    output.write(out_text)

# Save matrix as npz
save_npz('all_data_matrix.mtx.npz', all_data_subset, compressed=True)

# Save matrix as 10X Genomics mtx file
analysis.subset.csr_to_mtx(all_data_subset, out_name='all_data_matrix.mtx.gz', max_cells=100000)
```

Latly, if we want to reuse the PCA embeddings corrected during the integration step, we'll need to subset these as well. First, we'll load and subset the PCA data from the original 'analysis' object

```
_, pca_data, _, _, _ = analysis.load_reduced_dimensions(analysis.lower_dimensions_dir)
pca_data_subset = pca_data[cell_filter,].copy()
```

Then we'll add the pca_data_subset array to the analysis_subset object and invoke the reduce_dimensions function for Umap embedding and saving the reduced dimensions data and models locally.

```
analysis_subset.pca_data = pca_data_subset.copy()
analysis_subset.reduce_dimensions(use_integrated_pca=True, integrated_pca_path='', neighbors=30, n_permutations=3)
```

## 4 - Data (re)loading

The pipeline is modular, with different pieces of data saved to individual directories. This makes it easy to reload only what is needed or modify individual components.
The scripts below show an example of how to reload individual components (note that file names are the default).

```
# Init class
analysis = isa.integrated_analysis(max_cells_in_memory=200000)

# Load a merged dataset
analysis.load_merged_datasets('/path/to/4_merged_data')

# Set path to directory containing lower dimensions embeddings and models 
analysis.lower_dimensions_dir = '/path/to/5_lower_dimensions'

# Load HVGs list (saved in a text file, one gene per line)
analysis.load_hvgs('/path/to/common_higly_variable_genes.tsv.gz')

# Load clusters
analysis.load_clusters('/path/to/cluster_embeddings.tsv.gz')

# Load clusters markers
analysis.load_cluster_markers('/path/to/cluster_markers.tsv.gz')

# Set output directory for Umap plots
analysis.plot_dir = './'

# Load minimal spanning tree (i.e. trajectory)
analysis.load_trajectory_branches('/path/to/trajectory_analysis.txt')

# Load cells' pseudotime
analysis.load_trajectory_pseudotime('/path/to/trajectory_analysis_pseudotime.tsv.gz')

# Load cell cycle scores
analysis.load_cell_cycle_scores('/path/to/cell_cycle_scores.tsv.gz')
```
