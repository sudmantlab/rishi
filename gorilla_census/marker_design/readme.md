This repo explains how I produced the vcf file of putative target sites for carrying out targetted SNP sampling of mountain gorillas
Each markdown file explains a specific step in the process:
- [00_data_download.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/00_data_download.md) includes information on previously published data download (from SRA)
- [01_mapping.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/01_mapping.md) explains mapping these data and filtering the resulting bam files
- [02_genotyping.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/02_genotyping.md) covers genotyping and vcf file filtering (for quality etc.)
- [03_PCA.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/03_PCA.md) explains the production of PCAs from this data
- [04_geno.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/04_geno.md) outlines how to convert the vcf files to .geno files and then carry out popgen analyses on these files
- [05_mpcrselect](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/05_mpcrselect.md) covers running the mPCRselect pipeline to select putative target SNPs
- [06_snp_panel_popgen.md](https://github.com/sudmantlab/rishi/blob/main/gorilla_census/marker_design/06_snp_panel_popgen.md) plotting PCA from mpcrselect vcf output
