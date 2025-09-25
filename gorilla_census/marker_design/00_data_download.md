Code to download data to make gorilla SNP panel

move data to scratch dir (infiniband) using globus
to: `/global/scratch/users/rdekayne`  

make directory
```
cd /global/scratch/users/rdekayne/
mkdir -p /global/scratch/users/rdekayne/envs/
```

prepare conda environment for download
```
module load anaconda3
conda create -p /global/scratch/users/rdekayne/envs/sra
source activate base
conda activate /global/scratch/users/rdekayne/envs/sra
conda install bioconda::sra-tools
```

