/*
Workflow for for preprocessing scATAC-Seq data
*/

// ----------------Workflow---------------- //

include { BedtoolsMerge } from '../../modules/local/bedtools/bedtools_merge.nf'
include { AddSampleIDToBarcodes } from '../../modules/local/barcodes/add_sample_id_to_barcodes.nf'
include { ConcatenateBarcodes } from '../../modules/local/barcodes/concatenate_barcodes.nf'
include { BedtoolsCoverage } from '../../modules/local/bedtools/bedtools_coverage.nf'
include { AggregateCells } from '../../modules/local/aggregate/aggregate_cells.nf'

workflow SC_ATAC_PREPROCESSING {

  main:

  // CREATING INPUT CHANNELS -------------- //

  Channel
    .fromPath("${params.samples_manifest_file}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> file(row.peaks_bed) }
    .collect()
    .set{ bed_files }

  Channel
    .fromPath("${params.samples_manifest_file}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, file(row.barcodes)) }
    .set{ barcodes }

  Channel
    .fromPath("${params.samples_manifest_file}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, file(row.fragments)) }
    .set{ fragments }

  // CONSENSUS BED FILE ------------------- //

  BedtoolsMerge(bed_files)

  consensus_peaks = BedtoolsMerge.out.consensus_peaks

  // PROCESS BARCODES --------------------- //

  // Add sample ID information to barcode files

  AddSampleIDToBarcodes(barcodes)

  barcodes_with_sample_id = AddSampleIDToBarcodes.out.barcodes_with_sample_id.collect()

  // Concatenate the barcodes from all samples

  ConcatenateBarcodes(barcodes_with_sample_id)

  concatenated_barcodes = ConcatenateBarcodes.out.concatenated_barcodes

  // Convert the concatenated barcode file to a tuple channel

  concatenated_barcodes
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, row.barcode) }
    .set{ processed_barcodes }

  // COUNT PEAKS FRAGMENTS ---------------- //

  // Combine barcodes and fragments channel

  processed_barcodes
    .combine(fragments, by: 0)
    .set{ barcodes_and_fragments }

  // Compute peaks coverage for each cell

  BedtoolsCoverage(consensus_peaks, barcodes_and_fragments)

  // Group coverage files of cells from the same sample

  BedtoolsCoverage.out.peak_counts
    .groupTuple(by: 0)
    .set{ peak_counts }

  // AGGREGATE CELLS ---------------------- //

  // Channel for python script used for aggreagation

  aggregation_script = Channel.fromPath("${projectDir}/scripts/aggregate_peaks_counts.py")
  
  // Aggregate coverages of cells from the same sample, then save as MEX

  AggregateCells(aggregation_script, peak_counts)

}