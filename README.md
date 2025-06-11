# scRNA-Seq analysis tools
## Python modules for the analysis of single-cell RNA-Seq experiments

/// --------------------------------------- ///

### MODULES

#### 10X_analysis:

	Python module for the analysis of 10X scRNASeq experiments processed with CellRanger
	A tutorial is available

/// ---------- ///

#### integrated_analysis

	Python module for the integration and analysis of published datasets

### gene_clustering

	Python module for the clustering of genes in PCA space using Leiden algorithm
	(i.e. same procedure used to cluster cells in 10X_analysis and integrated_analysis modules)

/// --------------------------------------- ///

### NOTES

- 10X_analysis and integrated_analysis are designed for RNA Gene Expression data, but can handle ATAC data as well.

- For cell cycle analysis, the cell_cycle_genes_homo_sapiens.tsv and cell_cycle_genes_mus_musculus.tsv files report cell cycle gene sets for human and mouse data, respectively;

- For cluster identity analysis using the score_cluster_identity function, two tsv tables must be provided:

	- ref_samples table, with columns:

		- "short_name": short form name of the specific cell subclass (e.g. "B.Mem.Sp.v2");
		- "long_name": long form name of the specific cell subclass (e.g. "Spleen B cells");
		- "reference_cell_type": name of the main cell class (e.g. "B cells");

	- ref_expression table, with columns:

		- "GeneName": first column, reporting gene names;
		- Following columns are called as the elements from the "short_name" column in ref_samples table. Each row reports the expression value for the gene whose name is in the "GeneName" column;

- Import and pre-process files above as follows:

```
	# Import data
	ref_samples = pandas.read_csv('ref_samples.tsv', sep='\t')
	ref_expression = pandas.read_csv('ref_expression.tsv', sep='\t')

	# Calculate references Log2FC compared to all the other references
	gene_means = ref_expression.iloc[:, 1:].mean(axis=1)
	ref_expression.iloc[:, 1:] = numpy.log2(ref_expression.iloc[:, 1:].divide(gene_means, axis=0))
```