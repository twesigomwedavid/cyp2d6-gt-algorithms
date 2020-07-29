#!/usr/bin/env nextflow

bam_ch = Channel.fromFilePairs(params.in_pat) { file -> file.name.replaceAll(/.bam|.bai$/,'') }


process aldy {
  errorStrategy 'ignore'
     
  publishDir params.aldy_alleles_dir

  input:
     set val(name), file(bam) from bam_ch

  output:
     set val(name), file("${name}_2d6.aldy"), file("${name}_2d6.log") into aldy_ch

  script:
    
    """
        aldy genotype -p illumina -g cyp2d6 ${bam.get(0)} -o ${name}_2d6.aldy -l ${name}_2d6.log
    """

}
