Make PCAs from combined vcf

Make directory and conda environments
```
conda create -p /global/scratch/users/rdekayne/envs/pca
conda activate /global/scratch/users/rdekayne/envs/pca
conda install bioconda::plink

mkdir -p /global/scratch/users/rdekayne/gorilla_census/03_PCA && cd /global/scratch/users/rdekayne/gorilla_census/03_PCA
```
Now run plink to create PCAs which we plot later R - `03.1_PCA_unfilt.sh`
```
#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --time=0-4:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=pca.%j.out # output file
#SBATCH --error=pca.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=1 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/03_PCA

plink --vcf /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out 11gorillas_unfilt_out

touch /global/scratch/users/rdekayne/gorilla_census/03_PCA/11_pca.done
```
Now run - `sbatch 03.1_PCA_unfilt.sh`
We get: 3823659 variants and 14 people pass filters and QC.

And then download the data from the cluster `scp rdekayne@hpc.brc.berkeley.edu:/global/scratch/users/rdekayne/gorilla_census/03_PCA/*.eigen* .`
