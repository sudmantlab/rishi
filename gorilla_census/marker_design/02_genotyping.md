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

Make directory for filtered vcfs
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs
```
Then get a list of autosomes since we need to filt autosomes and sex chromosomes separately
```
cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | grep -v "chrX" | grep -v "chrY" > autosomes_scaffold_list.txt
```
Now we will filter the raw vcf files we made above - `02.5_vcf_filt2.sh`
Filter genotypes with low depth or low quality using `-e 'FORMAT/DP < 7 | FORMAT/GQ < 30'`
Set missing genotypes to ./. `--set-GTs .`
Exclude variant sites where the allele number (AN) is less than 2 `-e 'AN < 2'` i.e. filter out sites where all individuals are ./.
```
#!/bin/bash
#SBATCH --job-name=filt
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=filt.%j.out
#SBATCH --error=filt.%j.err

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/autosomes_scaffold_list.txt | sed -n ${numb}p)

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}"

bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz -e 'FORMAT/DP < 7 | FORMAT/GQ < 30' --set-GTs . -O u | bcftools filter -e 'AN < 2' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/"${scaf_name}".reheader_filt_mindepth7_minqual30.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/"${scaf_name}"_filt2.done
```
Submit `sbatch --array=1-23 02.5_vcf_filt2.sh`

Make list of vcf files to merge
```
ls /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/*.vcf.gz > autosomes_filt_to_merge.txt
```
Now concatenate all the autosomes into a single large vcf file - `02.6_vcf_auto_concat_filt.sh`
And then filter the resulting file to keep sites that are biallelic `-m2 -M2` have a minor allele count of 1 and at least 22 alleles `-e 'INFO/MAC < 1 | AN < 22'`
```
#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=concat.%j.out
#SBATCH --error=concat.%j.err

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge
bcftools concat -Ou -f /global/scratch/users/rdekayne/gorilla_census/02_genotyping/autosomes_filt_to_merge.txt | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30.vcf.gz

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2
#max 0 missing filter for SNPs
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30.vcf.gz -e 'INFO/MAC < 1 | AN < 22' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22.vcf.gz

touch concat_auto.done
```
Submit - `sbatch 02.6_vcf_auto_concat_filt.sh`

Count the number of sites retained
```
bcftools view -H autosomes_output_filt_mindepth7_minqual30_AN22.vcf.gz | wc -l
```
Results in: 4799902

Now do the same for sex chromosomes - `02.7_vcf_sex_filt.sh`
```
#!/bin/bash
#SBATCH --job-name=filt
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=filt.%j.out
#SBATCH --error=filt.%j.err

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtX
bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/chrX_mat_hsaX.reheader.raw.vcf.gz -e 'FORMAT/DP < 7 | FORMAT/GQ < 30' --set-GTs . -O u | bcftools filter -e 'AN < 2' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtX -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrX_mat_hsaX.reheader_filt_mindepth7_minqual30.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/chrX_mat_hsaX_filt.done

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtX2
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrX_mat_hsaX.reheader_filt_mindepth7_minqual30.vcf.gz -e 'INFO/MAC < 1 | AN < 19' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtX2 -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/chrX_mat_hsaX_output_filt_mindepth7_minqual30_AN19.vcf.gz


mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtY
bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/02_genotyping/raw_reheader_vcfs/chrY_pat_hsaY.reheader.raw.vcf.gz -e 'FORMAT/DP < 7 | FORMAT/GQ < 30' --set-GTs . -O u | bcftools filter -e 'AN < 2' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtY -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrY_pat_hsaY.reheader_filt_mindepth7_minqual30.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/02_genotyping/done_files/chrY_pat_hsaY_filt2.done

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtY2
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/02_genotyping/filt_vcfs/chrY_pat_hsaY.reheader_filt_mindepth7_minqual30.vcf.gz -e 'INFO/MAC < 1 | AN < 3' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_concat_filtY2 -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/chrY_pat_hsaY_output_filt_mindepth7_minqual30_AN06.vcf.gz
```
Run - `sbatch 02.7_vcf_sex_filt.sh`

Need to change the chromosome names to run mpcrselect (the target selection pipeline) `chr_name_conversion.txt`
```
chr1_pat_hsa1   chr1pathsa1
chr2_pat_hsa3   chr2pathsa3
chr3_pat_hsa4   chr3pathsa4
chr4_pat_hsa17x5    chr4pathsa17x5
chr5_mat_hsa6   chr5mathsa6
chr6_mat_hsa7   chr6mathsa7
chr7_pat_hsa8   chr7pathsa8
chr8_pat_hsa10  chr8pathsa10
chr9_pat_hsa11  chr9pathsa11
chr10_mat_hsa12 chr10mathsa12
chr11_mat_hsa2b chr11mathsa2b
chr12_pat_hsa2a chr12pathsa2a
chr13_pat_hsa9  chr13pathsa9
chr14_pat_hsa13 chr14pathsa13
chr15_pat_hsa14 chr15pathsa14
chr16_pat_hsa15 chr16pathsa15
chr17_mat_hsa18 chr17mathsa18
chr18_pat_hsa16 chr18pathsa16
chr19_pat_hsa5x17   chr19pathsa5x17
chr20_mat_hsa19 chr20mathsa19
chr21_pat_hsa20 chr21pathsa20
chr22_mat_hsa21 chr22mathsa21
chr23_mat_hsa22 chr23mathsa22
```
Change the name `02.8_vcf_mpcrselect_prep.sh`
```
#!/bin/bash
#SBATCH --job-name=prep
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=prep.%j.out # output file
#SBATCH --error=prep.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

bcftools annotate --rename-chrs chr_name_conversion.txt /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22.vcf.gz | bcftools sort -Oz -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.vcf.gz

bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.vcf.gz -e 'FORMAT/DP < 7 | (FORMAT/GQ) < 30' | bcftools sort -Oz -o /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename_refiltGQ30.vcf.gz

bcftools view -H /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.vcf.gz | wc -l
bcftools view -H /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename_refiltGQ30.vcf.gz | wc -l

touch full_filt.done
```
And run - `sbatch 02.8_vcf_mpcrselect_prep.sh`
