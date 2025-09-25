Code to create .geno files from .vcf files to calculate popgen statistics using [genomics_general](https://github.com/simonhmartin/genomics_general) tools
This specific section was carried out on an earlier iteration of the pipeline where we had 14 samples and is only retained here for clarity (hence why samples have AN28 rather than AN22 which was in the final version)

Download the set of tools from github and make directories
```
git clone https://github.com/simonhmartin/genomics_general

mkdir -p /global/scratch/users/rdekayne/gorilla_census/04_geno && cd /global/scratch/users/rdekayne/gorilla_census/04_geno
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos
```

Convert raw vcf files to geno files by filtering for no missing data i.e. AN = 28 - `04.1_geno_convert.sh`
```
#!/bin/bash
#SBATCH --job-name=geno
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=geno.%j.out # output file
#SBATCH --error=geno.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

module load python/3.11.6-gcc-11.4.0 bio/bcftools/1.16-gcc-11.4.0

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/autosomes_scaffold_list.txt | sed -n ${numb}p)

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}"

bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz --set-GTs . -e 'AN < 2' | bcftools filter -e 'INFO/MAC<1 | AN<28' | bcftools view -m1 -M2 -v snps | bcftools filter -e '(FORMAT/DP)<7 | (FORMAT/GQ)<30' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/"${scaf_name}"_filt_mindepth7_minqual30_m1M2_AN28.vcf.gz 

python /global/scratch/users/rdekayne/gorilla_census/04_geno/genomics_general/VCF_processing/parseVCF.py -i /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/"${scaf_name}"_filt_mindepth7_minqual30_m1M2_AN28.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=7 max=100 -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/"${scaf_name}"_filt_mindepth7_minqual30_m1M2_AN28.geno.gz

touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/"${scaf_name}".done
```
And run: `sbatch --array=1-23 04.1_geno_convert.sh`

Do the same with sex chromosomes:
cat ../02_genotyping/X_ploidy.txt | sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/01_mapping\/indiv_bams\///g' | sed 's/_1file.bam//g' > x_geno_ploidy.txt
cat ../02_genotyping/Y_ploidy.txt | sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/01_mapping\/indiv_bams\///g' | sed 's/_1file.bam//g' > y_geno_ploidy.txt
#error with y so put ploidy to 1

`04.2_geno_convert_sex.sh`
```
#!/bin/bash
#SBATCH --job-name=geno
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=geno.%j.out # output file
#SBATCH --error=geno.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

module load python/3.11.6-gcc-11.4.0 bio/bcftools/1.16-gcc-11.4.0

python /global/scratch/users/rdekayne/gorilla_census/04_geno/genomics_general/VCF_processing/parseVCF.py -i /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrX_mat_hsaX_filt_mindepth7_minqual30_m2M2_AN22.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=7 max=100 --ploidyFile x_geno_ploidy.txt -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/chrX_mat_hsaX_filt_mindepth7_minqual30_m2M2_AN22.geno.gz

python /global/scratch/users/rdekayne/gorilla_census/04_geno/genomics_general/VCF_processing/parseVCF.py -i /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrY_pat_hsaY_mindepth7_minqual30_m2M2_AN6.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=7 max=100 --ploidyFile y_all1_geno_ploidy.txt -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/chrY_pat_hsaY_mindepth7_minqual30_m2M2_AN6.geno.gz

touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/sex.done
```
And run - `sbatch 04.2_geno_convert_sex.sh`

Make a population file for calculating pop gen metrics like pi etc. 
`gorilla_popsfile.txt`
```
ind01_Maisha  mount
ind02_Tuck  mount
ind03_Turimaso  mount
ind04_Umurimo  mount
ind05_Imfura  mount
ind06_Kaboko  mount
ind07_Zirikana  mount
ind08_Dunia  elow
ind09_Itebero  elow
ind10_Pinga  elow
ind12_Tumani  elow
ind13_Ntabwoba  elow
ind14_Mkubwa  elow
ind15_Kaisi  elow
```
Now prep list for calculating pi per chromosome
```
ls /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/*AN28.vcf.gz > filt_vcfs_to_convert.txt
sed -i 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/02_genotyping\/filt_vcfs\///g' filt_vcfs_to_convert.txt
sed -i 's/_filt_mindepth7_minqual30_m2M2_AN28.vcf.gz//g' filt_vcfs_to_convert.txt
```
Now calculate pi in a loop - `04.3_geno_pi.sh`
Worth noting that genomics.py needs to be here
```
#!/bin/bash
#SBATCH --job-name=geno
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=geno.%j.out # output file
#SBATCH --error=geno.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

module load python/3.11.6-gcc-11.4.0 bio/bcftools/1.16-gcc-11.4.0

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/04_geno/filt_vcfs_to_convert.txt | sed -n ${numb}p)

python /global/scratch/users/rdekayne/gorilla_census/04_geno/genomics_general/popgenWindows.py -w 50000 -m 2000 -g /global/scratch/users/rdekayne/gorilla_census/02_genotyping/genos/"${scaf_name}"_filt_mindepth7_minqual30_AN28.geno.gz -o "${scaf_name}".div.regions.output.csv.gz -f phased -T 1 -p mount -p elow --popsFile /global/scratch/users/rdekayne/gorilla_census/04_geno/gorilla_popsfile.txt

touch /global/scratch/users/rdekayne/gorilla_census/04_geno/"${scaf_name}".pi.regions.done
```
Submit script `sbatch --array=1-23 04.3_geno_pi.sh`

Make directory for output and download pi files
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/04_geno/pi_output
cp *.div.regions.output.csv.gz /global/scratch/users/rdekayne/gorilla_census/04_geno/pi_output

scp rdekayne@hpc.brc.berkeley.edu:/global/scratch/users/rdekayne/gorilla_census/04_geno/pi_output/*.csv.gz .
```
