Code for the genotyping of gorilla loci for target SNP selection

Make directories and prepare conda environments
```
conda create -p /global/scratch/users/rdekayne/envs/geno
conda activate /global/scratch/users/rdekayne/envs/geno
conda install bioconda::bcftools

#make genotyping directory
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/ && cd /global/scratch/users/rdekayne/gorilla_census/02_genotyping/
#and make a directory for the done-files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files
#and make directory for raw vcf files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_vcfs
```

Want to get reference genome scaffolds in ascending order:
```
cut -f1 /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta.gz.fai > /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list
```
And get a list of bam files
```
ls /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/*.bam > bam.list
```

Get rid of 08-16 for low coverage and being eastern lowland
```
cat bam.list | grep -v 'ind08' |  grep -v 'ind09' |  grep -v 'ind10' |  grep -v 'ind11' |  grep -v 'ind12' | grep -v 'ind13' |  grep -v 'ind14' |  grep -v 'ind15' |  grep -v 'ind16' > 11_bam_list.txt
```

Now call genotypes using bcftools mpileup `02.1_genotyping.sh`
```
#!/bin/bash
#SBATCH --job-name=geno
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --output=geno.%j.out # output file
#SBATCH --error=geno.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job
#SBATCH -w n0032.savio4

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat w | sed -n ${numb}p)

bcftools mpileup -O u -d 250 --skip-indels -f /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta --annotate FORMAT/DP --bam-list /global/scratch/users/rdekayne/gorilla_census/02_genotyping/11_bam_list.txt -r "${scaf_name}" | bcftools call -m -f GQ -O v | bgzip > /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/"${scaf_name}".raw.vcf.gz.done
```
And submit `sbatch --array=1-23%1 02.1_genotyping.sh`

To genotype sex chromosomes we can specify ploidy - `ploidy.txt`
```
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind01_Maisha_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind02_Tuck_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind03_Turimaso_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind04_Umurimo_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind05_Imfura_1file.bam  M
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind06_Kaboko_1file.bam  M
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind07_Zirikana_1file.bam  M
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind17_Bwiruka_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind19_Katungi_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind20_Nyamunwa_1file.bam  F
/global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/ind21_Semehe_1file.bam  F
```
Then combine into male and female files
```
cat ploidy.txt | sed 's/ F/ 2/g' | sed 's/ M/ 1/g' > X_ploidy.txt
cat ploidy.txt | sed 's/ F/0/g' | sed 's/ M/ 1/g' > Y_ploidy.txt
```

And genotype these in a similar way using ploidy - `02.2_genotyping.sh`
```
#!/bin/bash
#SBATCH --job-name=geno_s
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --output=geno_s.%j.out # output file
#SBATCH --error=geno_s.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job
#SBATCH -w n0032.savio4

bcftools mpileup -O u -d 250 --skip-indels -f /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta --annotate FORMAT/DP --bam-list /global/scratch/users/rdekayne/gorilla_census/02_genotyping/11_bam_list.txt -r chrX_mat_hsaX | bcftools call -m --samples-file X_ploidy.txt -f GQ -O v | bgzip > /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_vcfs/chrX_mat_hsaX.raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/chrX_mat_hsaX.raw.vcf.gz.done

bcftools mpileup -O u -d 250 --skip-indels -f /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta --annotate FORMAT/DP --bam-list /global/scratch/users/rdekayne/gorilla_census/02_genotyping/11_bam_list.txt -r chrY_pat_hsaY | bcftools call -m --samples-file Y_ploidy.txt -f GQ -O v | bgzip > /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_vcfs/chrY_pat_hsaY.raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/chrY_pat_hsaY.raw.vcf.gz.done
```
Submit with `sbatch 02.2_genotyping.sh`

This produces vcf files but with weird naming so we want to reaheader the vcf files with short names (removing the filepath)
Get the current names
```
bcftools query -l raw_vcfs/chr10_mat_hsa12.raw.vcf.gz  > old_names.txt
```
Now make a new name by removing the filepath and ending
```
sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/01_mapping\/indiv_bams\///g' old_names.txt > new_names.txt
sed -i 's/_1file.bam//g' new_names.txt
```
Combine into one file
```
paste old_names.txt new_names.txt > reheader.txt
```
Make a new directory for the reheadered vcf files
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs
```

Now reheader them - `02.3_genotyping_rehead.sh`
```
#!/bin/bash
#SBATCH --job-name=rehead
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --output=rehead.%j.out # output file
#SBATCH --error=rehead.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

bcftools reheader --samples /global/scratch/users/rdekayne/gorilla_census/02_genotyping/reheader.txt -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/"${scaf_name}".rehead.done
```
Submit with `sbatch --array=1-25 02.3_genotyping_rehead.sh`

Index the resulting vcf files - `02.4_vcf_index.sh`
```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --output=index.%j.out # output file
#SBATCH --error=index.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

bcftools index /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz

##run
sbatch --array=1-25 02.4_vcf_index.sh
```














