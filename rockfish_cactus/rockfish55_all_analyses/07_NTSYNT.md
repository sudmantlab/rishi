
mkdir -p /global/scratch/users/rdekayne/rockfish55/07_NTSYNT && cd /global/scratch/users/rdekayne/rockfish55/07_NTSYNT
conda create -p /global/scratch/users/rdekayne/envs/synt
conda activate /global/scratch/users/rdekayne/envs/synt
conda install bioconda::mash
conda install -c bioconda -c conda-forge ntsynt
conda install --yes -c conda-forge -c bioconda quicktree snakemake intervaltree r-base bioconductor-treeio r-ggpubr bioconductor-ggtree r-phytools r-dplyr r-argparse r-scales r-stringr
R -e 'install.packages(c("gggenomes"), repos = "https://cran.r-project.org")'
tar xvzf ntSynt-viz-1.0.0.tar.gz
export PATH=/global/scratch/users/rdekayne/rockfish55/07_NTSYNT/ntSynt-viz-1.0.0/bin/:$PATH

ls ../assemblies/*_1.fasta > hap1_fasta_list.txt

#run_ntsynt.sh
#!/bin/bash
#SBATCH --job-name=ntsynt
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=ntsynt%j.out # output file
#SBATCH --error=ntsynt%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

ntSynt --fastas_list hap1_fasta_list.txt -d 7 -b 100,000 -p ntsynt_output_b10k_d7 -t 24

##run
sbatch run_ntsynt.sh

ls *.fai > fai_path_list.txt

#run_ntsynt_viz.sh
#!/bin/bash
#SBATCH --job-name=ntsynt
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=ntsynt%j.out # output file
#SBATCH --error=ntsynt%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

#ntsynt_viz.py --blocks ntsynt_output_b10k_d7.synteny_blocks.tsv --fais fai_path_list.txt --tree cactus_output_rooted_rearranged.nwk --target-genome SebastesalascanusSEB-2_shortspinethornyhead.fa.masked.hap1

ntsynt_viz.py --blocks ../ntsynt_output_b10k_d7.synteny_blocks.tsv --fais fai_path_list.txt --normalize --ribbon_adjust 0.7 --scale 1e9 --tree hap1_rockfish50_tree_rooted_rearranged.nwk --height 53 --width 60

##run
sbatch run_ntsynt_viz.sh
