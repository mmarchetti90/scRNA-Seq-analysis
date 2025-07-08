#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline for preprocessing scATAC-Seq data where peak calling was done on a per sample basis
*/

// ----------------Workflow---------------- //

include { SC_ATAC_PREPROCESSING } from './workflows/local/sc_atac_preprocessing.nf'

workflow {

  // Run workflow
  SC_ATAC_PREPROCESSING()

}