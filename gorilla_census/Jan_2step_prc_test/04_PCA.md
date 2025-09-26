Code for plotting PCAs from target snps in test run
Make directories and load conda env
```
conda activate /global/scratch/users/rdekayne/envs/pca

mkdir -p /global/scratch/users/rdekayne/gorilla_census/11_PCA && cd /global/scratch/users/rdekayne/gorilla_census/11_PCA
```

Now run PCA script `11.1_PCA_unfilt.sh`
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

plink --vcf /global/scratch/users/rdekayne/gorilla_census/10_genotyping/targets156_specific_filtered_MAC0_LRG.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out LRG_PCA_OUT

plink --vcf /global/scratch/users/rdekayne/gorilla_census/10_genotyping/targets156_specific_filtered_MAC0_SML.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out SML_PCA_OUT

touch /global/scratch/users/rdekayne/gorilla_census/11_PCA/pca.done
```
Submit `sbatch 11.1_PCA_unfilt.sh`
Output:
```
105 variants and 18 people pass filters and QC.
105 variants and 18 people pass filters and QC.
```
Adnd download the data
```
scp rdekayne@hpc.brc.berkeley.edu:/global/scratch/users/rdekayne/gorilla_census/11_PCA/*.eigen* .
```

