process ConcatenateBarcodes {

  // Concatenates barcode files
  
  label 'local'

  input:
  path(barcodes)

  output:
  path "all_barcodes.tsv", emit: concatenated_barcodes

  """
  echo -e "sample_id\tbarcode" > tmp.tsv
  cat *_barcodes.tsv >> tmp.tsv
  mv tmp.tsv all_barcodes.tsv
  """

}
