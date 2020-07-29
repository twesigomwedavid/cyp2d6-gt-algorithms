#!/bin/bash -l                                                                                    
#SBATCH -J vcf_qc                                                                              
#SBATCH -e errors_%j.txt                                                                          
#SBATCH -c 4 
#SBATCH --mem=8096 


ref="/path/to/reference.fa"


cd /path/to/vcfs/directory

for i in eur afr amr eas; do

    # Do joint calling using gvcfs

    gatk GenomicsDBImport -V getrm_${i}_gvcfs.list --genomicsdb-workspace-path getrm_${i}_database --intervals 22:42422500-42626900

    gatk GenotypeGVCFs -R $ref -V gendb://getrm_${i}_database -G StandardAnnotation -O getrm_${i}.vcf


    # Do hard filtering as region is too small for VQSR

    gatk SelectVariants -R $ref -V getrm_${i}.vcf -select-type SNP -O getrm_${i}_snps.vcf

    gatk VariantFiltration -R $ref -V getrm_${i}_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "check_snp" -O getrm_${i}_filtered_snps.vcf 

    gatk SelectVariants -R $ref -V getrm_${i}.vcf -select-type INDEL -O getrm_${i}_indels.vcf

    gatk VariantFiltration -R $ref -V getrm_${i}_indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "check_indel" -O getrm_${i}_filtered_indels.vcf 

done
