De-Kayne: Code to download data to make gorilla SNP panel

Make directory
```
cd /global/scratch/users/rdekayne/
mkdir -p /global/scratch/users/rdekayne/envs/
mkdir -p /global/scratch/users/rdekayne/gorilla_census && cd /global/scratch/users/rdekayne/gorilla_census
mkdir -p /global/scratch/users/rdekayne/gorilla_census/data && cd /global/scratch/users/rdekayne/gorilla_census/data
```

Prepare conda environment for download
```
module load anaconda3
conda create -p /global/scratch/users/rdekayne/envs/sra
source activate base
conda activate /global/scratch/users/rdekayne/envs/sra
conda install bioconda::sra-tools
```

`01_download_gorilla_data_01.sh`
```
#!/bin/bash
#SBATCH --job-name=sra
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl.out # output file
#SBATCH --error=dl.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

mkdir -p ind14_Mkubwa && cd ind14_Mkubwa
fasterq-dump --split-files SRX243452
touch SRX243452.done
fasterq-dump --split-files SRX243453
touch SRX243453.done
```






