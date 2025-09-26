Code for genotyping re-sequencing data for the 156 target SNPs  
Make directories and load conda envs
```
conda activate /global/scratch/users/rdekayne/envs/geno

#make genotyping directory
mkdir -p /global/scratch/users/rdekayne/gorilla_census/10_genotyping && cd /global/scratch/users/rdekayne/gorilla_census/10_genotyping
#and make a directory for the done-files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/10_genotyping/done_files
#and make directory for raw vcf files
mkdir -p /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_vcfs
```
Make a list of bam files for analysis
```
ls /global/scratch/users/rdekayne/gorilla_census/09_mapping/processed_bams/*.bam > new_bam_list.txt
```
Now do genotyping `10.1_genotyping.sh`
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

bcftools mpileup -O u -d 250 --skip-indels -f /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta --annotate FORMAT/DP --bam-list new_bam_list.txt -r "${scaf_name}" | bcftools call -m -f GQ -O v | bgzip > /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/10_genotyping/done_files/"${scaf_name}".raw.vcf.gz.done

```
Submit `sbatch --array=1-23 10.1_genotyping.sh`

Commands for reheadering these vcfs to remove filepath from the individual names
```
bcftools query -l raw_vcfs/chr10_mat_hsa12.raw.vcf.gz  > old_names.txt
sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/09_mapping\/processed_bams\///g' old_names.txt > new_names.txt
sed -i 's/.sorted.dup.bam//g' new_names.txt

paste old_names.txt new_names.txt > reheader.txt
```
Make directory for new vcfs
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_reheader_vcfs
```
Now reheader them `10.2_genotyping_rehead.sh`
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

bcftools reheader --samples /global/scratch/users/rdekayne/gorilla_census/10_genotyping/reheader.txt -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_vcfs/"${scaf_name}".reheader.raw.vcf.gz /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_vcfs/"${scaf_name}".raw.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_vcfs/"${scaf_name}".rehead.done

```
And run - `sbatch --array=1-25 10.2_genotyping_rehead.sh`

Now index these vcf files `10.3_vcf_index.sh`
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

bcftools index /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz
```
And run `sbatch --array=1-25 10.3_vcf_index.sh`

Make a dir for filtered vcfs
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/10_genotyping/filt_vcfs
```
Now filter the vcf files `10.4_vcf_filt2.sh`
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
scaf_name=$(cat l | sed -n ${numb}p)

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}"

#bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz -e 'AN < 2' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/filt_vcfs/"${scaf_name}".reheader_filt_AN2.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/10_genotyping/done_files/"${scaf_name}"_filt2.done

bcftools filter -Ou /global/scratch/users/rdekayne/gorilla_census/10_genotyping/raw_reheader_vcfs/"${scaf_name}".reheader.raw.vcf.gz -e 'AN < 1' | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_temp_"${scaf_name}" -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/filt_vcfs/"${scaf_name}".reheader_filt_AN1.vcf.gz && touch /global/scratch/users/rdekayne/gorilla_census/10_genotyping/done_files/"${scaf_name}"_filt2.done
```
And run: `sbatch --array=1-23 10.4_vcf_filt2.sh`

Make lists of files to concatenate (one for each filtering AN (missing data)
```
ls /global/scratch/users/rdekayne/gorilla_census/10_genotyping/filt_vcfs/*.vcf.gz > autosomes_filt_to_merge.txt
ls /global/scratch/users/rdekayne/gorilla_census/10_genotyping/filt_vcfs/*AN1.vcf.gz > autosomes_filt_to_mergeAN1.txt
```
Concatentate into a single file `10.5_vcf_auto_concat_filt.sh`
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
#bcftools concat -Ou -f /global/scratch/users/rdekayne/gorilla_census/10_genotyping/autosomes_filt_to_merge.txt | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz

#bcftools concat -Ou -f /global/scratch/users/rdekayne/gorilla_census/10_genotyping/autosomes_filt_to_mergeAN1.txt | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz

mkdir -p /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2
#max 0 missing filter for SNPs
#bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz -e 'INFO/MAC < 1' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz

bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz -e 'INFO/MAC < 1' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz

touch concat_auto.done
```
Run `sbatch 10.5_vcf_auto_concat_filt.sh`

Count loci in each concatenated file
```
bcftools view -H autosomes_output_filt_AN2_MAC1.vcf.gz | wc -l
```
321
```
bcftools view -H autosomes_output_filt_AN1_MAC1.vcf.gz | wc -l
```
321

Now filter for the exact amplicons we were targetting  

First index the files
```
bcftools index ./concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz
bcftools index ./concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz
```
Adjust list so we can filter
```
sed -i 's/pat/_pat_/g' target_markers.txt
sed -i 's/mat/_mat_/g' target_markers.txt
```
Then use bcftools -R for regions to filter and count loci
```
bcftools view -R target_markers.txt ./concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz > targets156_filtered.vcf
bcftools view -H targets156_filtered.vcf | wc -l
```
299

And for the other file
```
bcftools view -R target_markers.txt ./concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz > targets156_filteredAN1.vcf
bcftools view -H targets156_filteredAN1.vcf | wc -l
```
299

All the same counts!

Now do the same but instead of amplicons just extract specific SNP positions we were targetting
```
grep "_L" specific_targets.txt > specific_targets_filt.txt
sed -i 's/_L//g' specific_targets_filt.txt
sed -i 's/_/\t/g' specific_targets_filt.txt
sed -i 's/pat/_pat_/g' specific_targets_filt.txt
sed -i 's/mat/_mat_/g' specific_targets_filt.txt

bcftools view -R specific_targets_filt.txt ./concat_vcfs/autosomes_output_filt_AN2_MAC1.vcf.gz > targets156_specific_filtered.vcf
bcftools view -H targets156_specific_filtered.vcf | wc -l
```
85
```
bcftools view -R specific_targets_filt.txt ./concat_vcfs/autosomes_output_filt_AN1_MAC1.vcf.gz > targets156_specific_filteredAN1.vcf
bcftools view -H targets156_specific_filteredAN1.vcf | wc -l
```
85

Now we want to calulate missingness to see how well we captured variation
```
conda activate /global/scratch/users/rdekayne/envs/mpcrselectenv

vcftools --vcf targets156_specific_filtered.vcf --missing-indv
vcftools --vcf targets156_specific_filtered.vcf --missing-site
```
Now we will try without a MAC minor allele count filter `10.6_vcf_auto_concat_filt_inc_invariant.sh`
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
#bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN2.vcf.gz -e 'INFO/MAC < 0' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN2_MAC0.vcf.gz

bcftools filter -Oz /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN1.vcf.gz -e 'INFO/MAC < 0' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /global/scratch/users/rdekayne/gorilla_census/genotype_concat_merge2 -o /global/scratch/users/rdekayne/gorilla_census/10_genotyping/concat_vcfs/autosomes_output_filt_AN1_MAC0.vcf.gz

touch concat_auto.done
```
And submit `10.6_vcf_auto_concat_filt_inc_invariant.sh`

Check the output
```
bcftools view -H autosomes_output_filt_AN2_MAC0.vcf.gz | wc -l
#(OTHER MAC1) 321
# 485

bcftools view -H autosomes_output_filt_AN1_MAC0.vcf.gz | wc -l
#(OTHER MAC1) 321
# 485

#now filter for the list we wanted:
bcftools index ./concat_vcfs/autosomes_output_filt_AN2_MAC0.vcf.gz
bcftools index ./concat_vcfs/autosomes_output_filt_AN1_MAC0.vcf.gz
```

Now extract our specific loci from this file
```
bcftools view -R specific_targets_filt.txt ./concat_vcfs/autosomes_output_filt_AN2_MAC0.vcf.gz > targets156_specific_filtered_MAC0.vcf
bcftools view -H targets156_specific_filtered_MAC0.vcf | wc -l
```
105
```
bcftools view -R specific_targets_filt.txt ./concat_vcfs/autosomes_output_filt_AN1_MAC0.vcf.gz > targets156_specific_filtered_AN1MAC0.vcf
bcftools view -H targets156_specific_filtered_AN1MAC0.vcf | wc -l
```
105

Manage to get more loci this way (likeley some had missing data)

#now split our samples by the two different extraction kits we tired
```
bcftools view -S lrg_vol_sample_list.txt -o targets156_specific_filtered_MAC0_LRG.vcf targets156_specific_filtered_MAC0.vcf
bcftools view -S sml_vol_sample_list.txt -o targets156_specific_filtered_MAC0_SML.vcf targets156_specific_filtered_MAC0.vcf
```
Calculate missingness across kits
```
conda activate /global/scratch/users/rdekayne/envs/mpcrselectenv

vcftools --vcf targets156_specific_filtered_MAC0_LRG.vcf --missing-indv --out LRG_
vcftools --vcf targets156_specific_filtered_MAC0_SML.vcf --missing-indv --out SML_
```

Now cat the results `cat *_.imiss`
```
INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS
1-LargerVolumeStrip-well1_S1	105	0	41	0.390476
2-LargerVolumeStrip-well1_S2	105	0	39	0.371429
3-LargerVolumeStrip-well1_S3	105	0	33	0.314286
4-LargerVolumeStrip-well2_S4	105	0	86	0.819048
5-LargerVolumeStrip-well2_S5	105	0	101	0.961905
6-LargerVolumeStrip-well2_S6	105	0	86	0.819048
7-LargerVolumeStrip-well3_S7	105	0	61	0.580952
8-LargerVolumeStrip-well3_S8	105	0	66	0.628571
9-LargerVolumeStrip-well3_S9	105	0	72	0.685714
10-LargerVolumeStrip-well4_S10	105	0	45	0.428571
11-LargerVolumeStrip-well4_S11	105	0	58	0.552381
12-LargerVolumeStrip-well4_S12	105	0	70	0.666667
13-LargerVolumeStrip-well5_S13	105	0	61	0.580952
14-LargerVolumeStrip-well5_S14	105	0	61	0.580952
15-LargerVolumeStrip-well5_S15	105	0	100	0.952381
16-LargerVolumeStrip-well6_S16	105	0	103	0.980952
17-LargerVolumeStrip-well6_S17	105	0	103	0.980952
18-LargerVolumeStrip-well6_S18	105	0	101	0.961905
INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS
19-SmallerVolumeStrip-well1_S19	105	0	15	0.142857
20-SmallerVolumeStrip-well1_S20	105	0	15	0.142857
21-SmallerVolumeStrip-well1_S21	105	0	14	0.133333
22-SmallerVolumeStrip-well2_S22	105	0	54	0.514286
23-SmallerVolumeStrip-well2_S23	105	0	29	0.27619
24-SmallerVolumeStrip-well2_S24	105	0	43	0.409524
25-SmallerVolumeStrip-well3_S25	105	0	61	0.580952
26-SmallerVolumeStrip-well3_S26	105	0	60	0.571429
27-SmallerVolumeStrip-well3_S27	105	0	54	0.514286
28-SmallerVolumeStrip-well4_S28	105	0	55	0.52381
29-SmallerVolumeStrip-well4_S29	105	0	39	0.371429
30-SmallerVolumeStrip-well4_S30	105	0	57	0.542857
31-SmallerVolumeStrip-well5_S31	105	0	60	0.571429
32-SmallerVolumeStrip-well5_S32	105	0	29	0.27619
33-SmallerVolumeStrip-well5_S33	105	0	17	0.161905
34-SmallerVolumeStrip-well6_S34	105	0	96	0.914286
35-SmallerVolumeStrip-well6_S35	105	0	101	0.961905
36-SmallerVolumeStrip-well6_S36	105	0	100	0.952381
```
And in addition to individual missingness calculate site missingness 
```
vcftools --vcf targets156_specific_filtered_MAC0_LRG.vcf --missing-site --out LRG_
vcftools --vcf targets156_specific_filtered_MAC0_SML.vcf --missing-site --out SML_
```
LRG_.lmiss

SML_.lmiss










