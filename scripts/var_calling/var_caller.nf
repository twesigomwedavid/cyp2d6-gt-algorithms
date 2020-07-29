#!/usr/bin/env nextflow


bam_ch = Channel.fromFilePairs(params.in_align) { file -> file.name.replaceAll(/.bam|.bai$/,'') }

hg19ref = file(params.hg_ref, type: 'file')


process call_variants {

//   publishDir params.var_dir, mode: 'copy', overwrite: 'true'

//   label 'big_mem'
   input:
      set val(name), file(bam) from bam_ch

   output:
      set val(name), file("${name}.g.vcf") into gVCF_ch

   script:

    """
    gatk HaplotypeCaller -R ${hg19ref} -I ${bam.get(0)} -L 22:42422500-42626900 -ERC GVCF -stand-call-conf 30 -O ${name}.g.vcf

    """

}

process gzip {

   publishDir params.var_dir, mode: 'copy', overwrite: 'true'
   input:
      set val(name), file(gvcf) from gVCF_ch
   output:
      set val(name), file("${name}.g.vcf.gz"), file("${name}.g.vcf.gz.tbi") into final_output

   script:

    """
    bgzip ${gvcf}
    tabix ${name}.g.vcf.gz
    """

}