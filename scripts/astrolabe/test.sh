#!/bin/bash -l
#SBATCH -J astrolabe
#SBATCH -e errors_%j.txt
#SBATCH -c 4
#SBATCH --mem=20096


for i in $(cat sample_names.txt); do
    run-astrolabe.sh -inputVCF /path/to/vcfs/${i}.g.vcf.gz -skipBamQC -conf /opt/exp_soft/bioinf/astrolabe-0.8.7.0/astrolabe.ini -outFile /path/to/outfiles/${i}.log -verboseFile /path/to/alleles_dir/${i}.alleles.txt -novelFile /path/to/novel_dir/${i}.novel.txt -targets CYP2D6
done

