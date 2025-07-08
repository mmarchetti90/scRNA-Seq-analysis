process AddSampleIDToBarcodes {

  // Parses a barcode file and converts it to a tsv file with added sample info name
  
  label 'local'

  input:
  tuple val(sample_id), path(barcodes)

  output:
  path "${sample_id}_tagged_barcodes.tsv", emit: barcodes_with_sample_id

  """
  if [[ "${barcodes}" == *".gz" ]]
  then

    zcat ${barcodes} | awk -v FS='\t' -v OFS='\t' -v sample_id=${sample_id} '{ print sample_id, \$1 }' > ${sample_id}_tagged_barcodes.tsv

  else

    awk -v FS='\t' -v OFS='\t' -v sample_id=${sample_id} '{ print sample_id, \$1 }' ${barcodes} > ${sample_id}_tagged_barcodes.tsv

  fi
  """

}
