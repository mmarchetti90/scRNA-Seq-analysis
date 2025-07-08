process BedtoolsCoverage {

  // Get coverage of a cell's fragments on the consensus peaks
  
  label 'bedtools'

  input:
  each path(consensus_peaks)
  tuple val(sample_id), val(barcode), path(fragmets)

  output:
  tuple val(sample_id), path("${barcode}.cov.gz"), emit: peak_counts

  """
  # CellRanger does not count fragments overlapping with peaks
  # Instead, it counts the fragments ends (left and right) that fall within peaks
  # This because each end represents a transposition site, i.e. a chromatin accessible site

  if [[ "${fragmets}" == *".gz" ]]
  then

    zcat ${fragmets} | grep ${barcode} > tmp.fragments.bed

  else

    grep ${barcode} ${fragmets} > tmp.fragments.bed

  fi

  touch tmp.ends.bed
  awk -v FS='\t' -v OFS='\t' '{ print \$1,\$2,\$2+1 }' tmp.fragments.bed >> tmp.ends.bed
  awk -v FS='\t' -v OFS='\t' '{ print \$1,\$3,\$3+1 }' tmp.fragments.bed >> tmp.ends.bed

  sort -k1,1 -k2,2n tmp.ends.bed | gzip > ${barcode}.bed.gz

  bedtools coverage -a ${consensus_peaks} -b ${barcode}.bed.gz -counts | gzip > ${barcode}.cov.gz
  """

}
