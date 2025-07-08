#!/bin/bash

barcodes=${1}
peaks_fragment=${2}
consensus_peaks=${3}
sample_name="${4:-atac_data}"

### Init output dir

rm -r ${sample_name}
mkdir -p ${sample_name}/cell_frags
mkdir -p ${sample_name}/cell_covs

### Unzip barcodes

if [[ "${barcodes}" == *".gz" ]]
then

	cp ${barcodes} barcodes_copy.txt.gz
	gzip -d barcodes_copy.txt.gz
	barcodes=barcodes_copy.txt

fi

### Extract fragments and compute coverage for each cell

while IFS= read -r cell_barcode
do

	if [[ "${peaks_fragment}" == *".gz" ]]
	then

		zcat ${peaks_fragment} | grep ${cell_barcode} | gzip > ${sample_name}/cell_frags/${cell_barcode}.bed.gz

	else

		grep ${cell_barcode} ${peaks_fragment} | gzip > ${sample_name}/cell_frags/${cell_barcode}.bed.gz

	fi

	bedtools coverage -a ${consensus_peaks} -b ${sample_name}/cell_frags/${cell_barcode}.bed.gz -counts | gzip > ${sample_name}/cell_covs/${cell_barcode}.cov.gz

done < "${barcodes}"

rm barcodes_copy.txt
