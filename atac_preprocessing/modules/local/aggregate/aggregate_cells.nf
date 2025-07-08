process AggregateCells {

  // Aggregate ATAC peak counts from individual cells into a MEX matrix
  
  label 'python'

  publishDir "${params.main_output_dir}", mode: "copy", pattern: "${sample_id}"
  publishDir "${params.main_output_dir}", mode: "copy", pattern: "*_mex.stats"

  input:
  each path(aggregation_script)
  tuple val(sample_id), path(counts)

  output:
  tuple val(sample_id), path("${sample_id}"), emit: mex
  tuple val(sample_id), path("${sample_id}_mex.stats"), emit: mex_stats

  """
  python ${aggregation_script} \
  --sample_name ${sample_id} \
  --cov_dir . \
  --min_detected_fragments ${params.min_detected_fragments} \
  --min_detected_peaks ${params.min_detected_peaks} &> ${sample_id}_mex.stats
  """

}
