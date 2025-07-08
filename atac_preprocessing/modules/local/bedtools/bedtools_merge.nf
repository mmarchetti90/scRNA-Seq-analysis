process BedtoolsMerge {

  // Merge bed files to create a consensus
  
  label 'bedtools'

  publishDir "${params.main_output_dir}/${params.consensus_peaks_dir}", mode: "copy", pattern: "consensus_peaks.bed"

  input:
  path(bed_files)

  output:
  path "consensus_peaks.bed", emit: consensus_peaks

  """
  cat *.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > consensus_peaks.bed
  """

}
