# Comparison of Pharmacogene Star Allele Calling Algorithms

This repository contains the code used in the analysis for our [manuscript](https://www.nature.com/articles/s41525-020-0135-2):

Twesigomwe, D., Wright, G.E.B., Drögemöller, B.I. et al. A systematic comparison of pharmacogene star allele calling bioinformatics algorithms: a focus on CYP2D6 genotyping. npj Genom. Med. 5, 30 (2020). https://doi.org/10.1038/s41525-020-0135-2


## Abstract

Genetic variation in genes encoding cytochrome P450 enzymes has important clinical implications for drug metabolism. Bioinformatics algorithms for genotyping these highly polymorphic genes using high-throughput sequence data and automating phenotype prediction have recently been developed. The CYP2D6 gene is often used as a model during the validation of these algorithms due to its clinical importance, high polymorphism, and structural variations. However, the validation process is often limited to common star alleles due to scarcity of reference datasets. In addition, there has been no comprehensive benchmark of these algorithms to date. We performed a systematic comparison of three star allele calling algorithms using 4618 simulations as well as 75 whole-genome sequence samples from the GeT-RM project. Overall, we found that Aldy and Astrolabe are better suited to call both common and rare diplotypes compared to Stargazer, which is affected by population structure. Aldy was the best performing algorithm in calling CYP2D6 structural variants followed by Stargazer, whereas Astrolabe had limitations especially in calling hybrid rearrangements. We found that ensemble genotyping, characterised by taking a consensus of genotypes called by all three algorithms, has higher haplotype concordance but it is prone to ambiguities whenever complete discrepancies between the tools arise. Further, we evaluated the effects of sequencing coverage and indel misalignment on genotyping accuracy. Our account of the strengths and limitations of these algorithms is extremely important to clinicians and researchers in the pharmacogenomics and precision medicine communities looking to haplotype CYP2D6 and other pharmacogenes using high-throughput sequencing data.


*This code may need to be modified to work on your computational platform*
