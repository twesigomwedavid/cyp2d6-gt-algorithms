#!/usr/bin/env nextflow

// This script automates simulation of CYP2D6/2D7/2D8 region NGS data for random combinations of CYP2D6 SNV-defined star alleles. 

// Authors: David Twesigomwe and Phelelani Mpangase

// Date last modified: December 06, 2019


Channel.fromPath("/path/to/pharmvar/haplotype_vcfs/CYP2D6/GRCh37/*vcf.gz")
    .randomSample(160)
    .buffer (size: 2, remainder: true)
    //.println { "${it}" }
    .set { data }


data.into { data1; data2 }

runfile = file(params.run_file, type: 'file')
hg19ref = file(params.hg19, type: 'file')

//data2.subscribe { println "${it}" } 

process make_haplo {
    input:
    each pair from data1
    
    output:
    set val("sample"), file("*_{maternal,paternal}.fasta") into out_file
    
    """
    gatk FastaAlternateReferenceMaker -R ${hg19ref} -L 22:42518000-42555000 -V ${pair.get(0)} -O ${pair.get(0).name.replace(/.vcf.gz/, "")}_maternal.fasta
    gatk FastaAlternateReferenceMaker -R ${hg19ref} -L 22:42518000-42555000 -V ${pair.get(1)} -O ${pair.get(1).name.replace(/.vcf.gz/, "")}_paternal.fasta
    """
}

process sim_seq {
    input:
    set sample, file(mutated) from out_file
    
    output:
    set sample, file("*.fastq") into simulated

    """
    simLibrary -i 400 -x 8 ${mutated.get(0)} | simNGS -o fastq -p paired ${runfile} > ${mutated.find { it =~ '_maternal.fasta$' }.baseName }.simulated
    simLibrary -i 400 -x 8 ${mutated.get(1)} | simNGS -o fastq -p paired ${runfile} > ${mutated.find { it =~ '_paternal.fasta$' }.baseName }.simulated
    cat ${mutated.find { it =~ '_maternal.fasta$' }.baseName }.simulated ${mutated.find { it =~ '_paternal.fasta$' }.baseName }.simulated > ${mutated.find { it =~ '_paternal.fasta$' }.baseName }_${mutated.find { it =~ '_maternal.fasta$' }.baseName }.fastq
    """
}

//simulated.subscribe { println "${it}" }

process align {
//   maxForks 10
//   errorStrategy 'ignore'
//     tag "${sim_reads.baseName}"
   label 'big_mem'
   input:
   set sample, file(sim_reads) from simulated
   output:
   set sample, file("*.sam") into alignsam_ch
//   name = sim_reads.name.replace(/.fastq/, "")

   script: 

   """
   bwa mem -M -t 4 -p ${hg19ref} ${sim_reads} > ${sim_reads.name.replace(/.fastq/, "")}.sam
   """
}

//alignsam_ch.subscribe { println "${it}" }

process get_bam {
 //  maxForks 10

   input:
   set sample, file(sam) from alignsam_ch
   output:
   set sample, file("*_conv.bam") into alignbam_ch

   script:
     """
         samtools view -bS ${sam} -o ${sam.name.replace(/.sam/, "")}_conv.bam
     """

}

process add_read_groups {

   input:
      set sample, file(bam) from alignbam_ch

   output:
      set sample, file("*_addRG.bam") into newbam_ch

   script:

    """
    gatk AddOrReplaceReadGroups -I ${bam} -O ${bam.name.replace(/_conv.bam/, "")}_addRG.bam -ID SIM3.1 -LB library1 -PL illumina -PU SIM3RX30.1 -SM ${bam.name.replace(/_conv.bam/, "")} -SO coordinate
    """
}

process rename_bams {

   publishDir params.aln_dir, mode: 'copy', overwrite: 'true'

   input:
      set sample, file(bamRG) from newbam_ch

   output:
      set sample, file("*.{bam,bai}") into fixedbams_ch

   script:


    """
    mv ${bamRG} ${bamRG.name.replace(/_addRG.bam/, "")}.bam
    samtools index ${bamRG.name.replace(/_addRG.bam/, "")}.bam

    """

}
