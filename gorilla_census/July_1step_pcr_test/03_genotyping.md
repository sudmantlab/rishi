Code for genotyping of 1step PCR test

Load conda env and make directories
```
conda activate /global/scratch/users/rdekayne/envs/geno

#make genotyping directory
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping
#and make a directory for the done-files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/done_files
#and make directory for raw vcf files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_vcfs
```
List of bam files we will use is here: `/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/bam_96_list.txt`

Now do genotyping with `genotyping.sh`
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

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

bcftools mpileup -O u -d 250 --skip-indels -f /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta --annotate FORMAT/DP --bam-list /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/bam_96_list.txt -r "${scaf_name}" | bcftools call -m -f GQ -O v | bgzip > /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/done_files/"${scaf_name}".raw.vcf.gz.done
```
and run `sbatch --array=1-23 genotyping.sh`

We are going to reheader the vcfs so remove the filepath and prepare a file for the renaming
```
bcftools query -l raw_vcfs/chr10_mat_hsa12.raw.vcf.gz > old_names.txt
sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/July_2025_library_test_UCB\/processed_bams\///g' old_names.txt > new_names.txt
sed -i 's/.sorted.dup.bam//g' new_names.txt

paste old_names.txt new_names.txt > reheader.txt

mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_reheader_vcfs
```
Run VCF reheadering script `genotyping_rehead.sh`
```
#!/bin/bash
#SBATCH --job-name=rehead
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=rehead.%j.out
#SBATCH --error=rehead.%j.err

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

bcftools reheader --samples /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/reheader.txt -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_vcfs/"${scaf_name}".reheader.raw.vcf.gz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_vcfs/"${scaf_name}".rehead.done
```
Submit with `sbatch --array=1-23 genotyping_rehead.sh`

And move the new vcf files to the correct directory
```
mv raw_vcfs/*.reheader* raw_reheader_vcfs/
```

Index the vcf files `vcf_index.sh`
```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=index.%j.out
#SBATCH --error=index.%j.err

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

bcftools index /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz
```
Submit `sbatch --array=1-23 vcf_index.sh`

Now make a new directory for our filtered files
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs
```

Filter the VCF files `vcf_filt2.sh`
```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=index.%j.out
#SBATCH --error=index.%j.err

numb=${SLURM_ARRAY_TASK_ID}
scaf_name=$(cat /global/scratch/users/rdekayne/gorilla_census/02_genotyping/scaf.list | sed -n ${numb}p)

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}"

#bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz -e 'AN < 2' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/"${scaf_name}".reheader_filt_AN2.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/"${scaf_name}"_filt2.done

bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz -e 'AN < 1' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/"${scaf_name}".reheader_filt_AN1.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/"${scaf_name}"_filt22.done
```
And run `sbatch --array=1-23 vcf_filt2.sh`

Make list of files to merge from each filter
```
ls /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/*AN2.vcf.gz > autosomes_filt_to_mergeAN2.txt
ls /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/filt_vcfs/*AN1.vcf.gz > autosomes_filt_to_mergeAN1.txt
```

Concatenate vcfs into a single file for each filtering type with `vcf_auto_concat_filt.sh`
```
#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=concat.%j.out
#SBATCH --error=concat.%j.err

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge
bcftools concat -Ou -f /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/autosomes_filt_to_mergeAN2.txt | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz
mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge
bcftools concat -Ou -f /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/autosomes_filt_to_mergeAN1.txt | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2
#max 0 missing filter for SNPs
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz -e 'INFO/MAC < 1' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz
mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz -e 'INFO/MAC < 1' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz

touch concat_auto.done
```
And run `sbatch vcf_auto_concat_filt.sh`

Check the number of SNPs output
```
bcftools view -H concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz | wc -l
```
278
```
bcftools view -H concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz | wc -l
```
278

Also make file with no MAC filter to retain weird sites `vcf_auto_concat_filt_inc_invariant.sh`
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

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2
#max 0 missing filter for SNPs
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz -e 'INFO/MAC < 0' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o ./concat_vcfs/autosomes_output_filt_AN2_MAC0.vcf.gz
mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge3
bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz -e 'INFO/MAC < 0' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge3 -o ./concat_vcfs/autosomes_output_filt_AN1_MAC0.vcf.gz

touch concat_auto.done
```
And run `sbatch vcf_auto_concat_filt_inc_invariant.sh`

Now we want to filter vcfs to get only our 40 target sites

We have a file `specific_targets_filt.txt` containing these targets:
``` 
chr2_pat_hsa3	63820727
chr4_pat_hsa17x5	139673146
chr15_pat_hsa14	36212595
chr1_pat_hsa1	10233235
chr6_mat_hsa7	11390729
chr7_pat_hsa8	5143130
chr19_pat_hsa5x17	24137428
chr1_pat_hsa1	51992760
chr2_pat_hsa3	104430425
chr3_pat_hsa4	121313448
chr4_pat_hsa17x5	181243267
chr5_mat_hsa6	61734828
chr6_mat_hsa7	163795044
chr7_pat_hsa8	39191000
chr7_pat_hsa8	45067101
chr7_pat_hsa8	63184015
chr8_pat_hsa10	86733664
chr8_pat_hsa10	145958666
chr9_pat_hsa11	44720414
chr9_pat_hsa11	70173757
chr10_mat_hsa12	20613085
chr13_pat_hsa9	35940167
chr17_mat_hsa18	82026819
chr18_pat_hsa16	103784314
chr1_pat_hsa1	180923230
chr4_pat_hsa17x5	169976139
chr8_pat_hsa10	26547585
chr9_pat_hsa11	134720610
chr10_mat_hsa12	107788170
chr18_pat_hsa16	74778095
chr4_pat_hsa17x5	49217303
chr5_mat_hsa6	93647643
chr9_pat_hsa11	27699017
chr10_mat_hsa12	64331685
chr10_mat_hsa12	75066868
chr16_pat_hsa15	55694687
chr16_pat_hsa15	79107645
chr18_pat_hsa16	111163844
chr19_pat_hsa5x17	18421770
chr22_mat_hsa21	25719657
```

To do this we will use our unfiltered concatenated file
```
bcftools view -R specific_targets_filt.txt /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz > targets40_specific_filtered_AN2.vcf
bcftools view -H targets40_specific_filtered_AN2.vcf | wc -l
```
40
```
bcftools view -R specific_targets_filt.txt /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz > targets40_specific_filtered_AN1.vcf
bcftools view -H targets40_specific_filtered_AN1.vcf | wc -l
```
40

These both contain the same information and info of all 40 target snps so we will now plot these using an R script to show a heatmap

