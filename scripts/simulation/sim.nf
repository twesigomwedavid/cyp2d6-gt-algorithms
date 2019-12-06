#!/usr/bin/env nextflow

// This script generates simulated BAM files from haplotype sequences. It also adds alignments for control regions required for CYP2D6 SV calling using Astrolabe, Aldy, and Astrolabe.

// Author: David Twesigomwe

// Date last modified: December 6, 2019

// Input and output file parameters are set in the nextflow.config file 

input_ch = Channel.fromFilePairs(params.in_pat, type: 'file')

ctrl_gene_reads = file(params.ctrl_path, type: 'file')
astro_ctrl = file(params.astro_control, type: 'file')
runfile = file(params.run_file, type: 'file')
hg19ref = file(params.hg19, type: 'file')


process simulation {
//   maxForks 10
//   publishDir params.fastq_dir, mode: 'copy', overwrite: 'true'

   input:
      set val(name), file(fasta) from input_ch
   
   output:
      set val(name), file("${name}_final.fastq") into reads_ch
      
   script:
     """
         simLibrary -i 400 -x 8 ${fasta.get(0)} | simNGS -o fastq -p paired ${runfile} > ${name}_maternal.fastq
         simLibrary -i 400 -x 8 ${fasta.get(1)} | simNGS -o fastq -p paired ${runfile} > ${name}_paternal.fastq
	 simLibrary -i 400 -x 8 ${astro_ctrl} | simNGS -o fastq -p paired ${runfile} > astro_ct_maternal.fastq
	 simLibrary -i 400 -x 8 ${astro_ctrl} | simNGS -o fastq -p paired ${runfile} > astro_ct_paternal.fastq
         cat ${name}_maternal.fastq ${name}_paternal.fastq ${ctrl_gene_reads} astro_ct_maternal.fastq astro_ct_paternal.fastq > ${name}_final.fastq 
     """

}

process align {
//   maxForks 10
   errorStrategy 'ignore'
   tag "${name}"
   label 'big_mem'
   input:
      set val(name), file(fastq) from reads_ch
   output:
      set val(name), file("${name}.sam") into alignsam_ch
      
   script:
     """
     bwa mem -M -t 4 -p ${hg19ref} ${fastq} > ${name}.sam
     """

}

process get_bam {
 //  maxForks 10

   input:
      set val(name), file(sam) from alignsam_ch
   output:
      set val(name), file("${name}.bam") into alignbam_ch

   script:
     """
         samtools view -bS ${sam} -o ${name}.bam
     """

}

//alignbam_ch.subscribe { println "${it}" }


process add_read_groups {

   input:
      set val(name), file(bam) from alignbam_ch

   output:
      set val(name), file("${name}_addRG.bam") into newbam_ch

   script:

    """
    gatk AddOrReplaceReadGroups -I ${bam} -O ${name}_addRG.bam -ID SIM3.1 -LB library1 -PL illumina -PU SIM3RX30.1 -SM $name -SO coordinate
    """
}

//newbam_ch.subscribe {println "$it"}

process rename_bams {

   publishDir params.aln_dir, mode: 'copy', overwrite: 'true'

   input:
      set val(name), file(bamRG) from newbam_ch

   output:
      set val(name), file("${name}.bam"), file("${name}.bam.bai") into fixedbams_ch

   script:

    """
    mv ${bamRG} ${name}.bam
    samtools index ${name}.bam 

    """

}

