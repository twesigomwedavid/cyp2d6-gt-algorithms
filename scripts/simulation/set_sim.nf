#!/usr/bin/env nextflow

// This script generates CYP2D6/2D7/2D8 regional NGS data simulations for combinations of various CYP2D6 SNV-defined star alleles with a fixed haplotype. 

// Author: David Twesigomwe

// Date last modified: December 6, 2019

pgvcfs_ch = Channel.fromFilePairs(params.in_pgvcf, type: 'file') { file -> file.name.replaceAll(/.vcf.gz|.tbi$/,'') }

runfile = file(params.run_file, type: 'file')
hg19ref = file(params.hg19, type: 'file')
fixed_haplo = file(params.test_haplo, type: 'file')


process make_haplo {
   input:
	
      set val(name), file(pharmgvcf) from pgvcfs_ch
   output:	   
      set val(name), file("${name}_P.fasta") into seq_ch
      
   script:
    """ 
         gatk FastaAlternateReferenceMaker -R ${hg19ref} -L 22:42518000-42555000 -V ${pharmgvcf.get(0)} -O ${name}_P.fasta
    """

}

process simulation {
//   maxForks 10
   publishDir params.fastq_dir, mode: 'copy', overwrite: 'true'

   input:
      set val(name), file(fasta) from seq_ch
        
   output:
      set val(name), file("${name}_final.fastq") into reads_ch
      
   script:
     """ 
         simLibrary -i 400 -x 8 ${fasta} | simNGS -o fastq -p paired ${runfile} > ${name}_P.fastq
         cat ${fixed_haplo} ${name}_P.fastq > ${name}_final.fastq  
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


process rename_bams {

   publishDir params.aln_dir, mode: 'copy', overwrite: 'true'

   input:
      set val(name), file(bamRG) from newbam_ch

   output:
      set val(name), file("${name}_anystring.bam"), file("${name}_anystring.bam.bai") into fixedbams_ch

   script:

    """
    mv ${bamRG} ${name}_anystring.bam
    samtools index ${name}_anystring.bam 

    """

}

