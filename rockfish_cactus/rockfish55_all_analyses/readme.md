This folder contains all documents related to the analysis of 55 rockfish haplotypes done in 2025 by De-Kayne.

The steps outlined are as follows:

- [00_raw_data.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/00_raw_data.md) - this contains the code to compile the core assemblies that will be used in this study
- [01_BUSCOS.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/01_BUSCOS.md) - this contains code to run the BUSCO analysis identifying conserved genes across haplotype resolved assemblies
- [02_TREES.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/02_TREES.md) - here we produce trees from BUSCO genes identified preivously and produce a species tree that we will use for cactus input
- [03_CACTUS.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/03_CACTUS.md) - here is the code to run cactus and parse the output
- [04_HAL_TO_TREES.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/04_HAL_TO_TREES.md) - this code explains extracting 100kb window msas from our full .hal cactus output format and then producing and analysing trees from these regions
- [05_GENE_RER.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/05_GENE_RER.md) - next we produce trees from genes by extracting exon msas from our cactus output and then carry out RERconverge analysis on these gene trees
- [06_PSMC.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/06_PSMC.md) - here is the code to run psmc and plot the output on our different species
- [07_NTSYNT.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/07_NTSYNT.md) - here is the code to identify patterns of synteny across assemblies using NTSynt
- [08_ANNOTATIONS.md](https://github.com/sudmantlab/rishi/blob/main/rockfish_cactus/rockfish55_all_analyses/08_ANNOTATIONS.md) - here is code to parse and analyse CAT output from the annotation of all our haplotypes from the honeycomb genome annotation
