This repo explains how I produced the vcf file of putative target sites for carrying out targetted SNP sampling of mountain gorillas
Each markdown file explains a specific step in the process:
- 00_data_download.md includes information on previously published data download (from SRA)
- 01_mapping.md explains mapping these data and filtering the resulting bam files
- 02_genotyping.md covers genotyping and vcf file filtering (for quality etc.)
- 03_PCA.md explains the production of PCAs from this data
- 04_geno.md outlines how to convert the vcf files to .geno files and then carry out popgen analyses on these files
