In this section we will make trees from the busco genes we identified in 01

First make environments and prepare environment
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/02_TREES && cd /global/scratch/users/rdekayne/rockfish55/02_TREES
conda activate BUSCO_phylogenomics

#python BUSCO_phylogenomics/BUSCO_phylogenomics.py -help
cp -r /global/scratch/users/rdekayne/rockfish/01_BUSCOS/BUSCO_phylogenomics .
```

We will then run [busco_phylogenomics](https://github.com/jamiemcg/BUSCO_phylogenomics) with `run_BUSCO_phyhlogenomics.sh`
```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --output=busco_phylo%j.out # output file
#SBATCH --error=busco_phylo%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

python ./BUSCO_phylogenomics/BUSCO_phylogenomics.py  -i /global/scratch/users/rdekayne/rockfish55/01_BUSCOS/output_busco_phylogenomics/ -o ./complete_busco_phylogenomics --gene_trees_only --gene_tree_program iqtree -t 24

touch BUSCO_phylogenomics.done
```
and run this with `sbatch run_BUSCO_phyhlogenomics.sh`

Now prepare ASTRAL to run
```
conda install bioconda::astral-tree
git clone https://github.com/smirarab/ASTRAL/
module load openjdk/17.0.8.1_1-gcc-11.4.0
./make.sh
```

And make a directory for our output
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/02_TREES/astral_run && cd /global/scratch/users/rdekayne/rockfish55/02_TREES/astral_run
```

And then run the ASTRAL analysis `make_astral_tree.sh`
```
#!/bin/bash
#SBATCH --job-name=astral
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=astral%j.out # output file
#SBATCH --error=astral%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

java -jar /global/scratch/users/rdekayne/rockfish/01_BUSCOS/ASTRAL/astral.5.7.8.jar -i /global/scratch/users/rdekayne/rockfish55/02_TREES/complete_busco_phylogenomics/gene_trees_single_copy/ALL.tree -o all_buscos_55rockfish.tre
touch tree.done
```
and submit `sbatch make_astral_tree.sh`
