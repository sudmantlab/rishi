Code for mapping reads from 1step pcr test run

Make directory and load conda environment
```
mkdir /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/mapping && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/mapping
conda activate /global/scratch/users/rdekayne/envs/mapping

mkdir /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams
```
Make a list of sample prefixes
```
ls /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/combined_indiv_fastas | grep _R1.fastq.gz | sed 's/_R1.fastq.gz//g' > sample_prefixes.txt
```
Now run mapping script `mapping_map_p1.sh `
```
#!/bin/bash
#SBATCH --job-name=bwa_map
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=bwa_map.%j.out
#SBATCH --error=bwa_map.%j.err

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/mapping/sample_prefixes.txt | sed -n ${ind}p)

# call bwa
bwa-mem2 mem -t 24 /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/combined_indiv_fastas/${indiv_name}_R1.fastq.gz /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/combined_indiv_fastas/${indiv_name}_R2.fastq.gz | samtools sort -@24 -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/${indiv_name}.raw.bam && touch ${indiv_name}.mapping.done
```
And run `sbatch --array=1-96 mapping_map_p1.sh`

Make some directories, move error and out files and remove done files
```
mkdir all_err_out
mv *.err all_err_out/
mv *.out all_err_out/
rm *.done
```

Now for bam processing, load java and make directoires
```
module load java/22.0.1
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam
```
Now run the processing script `mapping_process_p1.sh `
```
#!/bin/bash
#SBATCH --job-name=bwa_proc
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=bwa_proc.%j.out
#SBATCH --error=bwa_proc.%j.err

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/mapping/sample_prefixes.txt | sed -n ${ind}p)

#PART1
picard FixMateInformation I=/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/${indiv_name}.raw.bam VALIDATION_STRINGENCY=LENIENT OUTPUT=/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam/${indiv_name}.raw.bam

sambamba sort /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam/${indiv_name}.raw.bam -o /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam/${indiv_name}.sorted.raw.bam -t 24 -m 50GB

picard MarkDuplicates INPUT=/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam/${indiv_name}.sorted.raw.bam OUTPUT=/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams/${indiv_name}.sorted.dup.bam METRICS_FILE=/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams/${indiv_name}.sorted.dup.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024

sambamba index /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams/${indiv_name}.sorted.dup.bam

#now remove all intermediate files
rm /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/raw_bams/tmp_bam/${indiv_name}*

touch ${indiv_name}.processing.done
```
Submit `sbatch --array=1-96%8 mapping_process_p1.sh  `

Move error and out files
```
mkdir all_err_out
mv *.err all_err_out/
mv *.out all_err_out/
```
And make sample list of bam files for calling genotypes
```
ls /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams/*.bam > bam_96_list.txt
```

Check the depth of specific target loci using mosdepth
```
conda activate /global/scratch/users/rdekayne/envs/mapping
```
Prepare bam list
```
sed 's/\/global\/scratch\/users\/rdekayne\/gorilla_census\/July_2025_library_test_UCB\/processed_bams\///g' bam_96_list.txt > bam_96_list_simplified.txt
```
Then prepare input bed file for mosdepth `mosdepth_40_targets.bed`
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/depth_calc && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/depth_calc
cp /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/genotyping/specific_targets_filt.txt .
awk -F'\t' -v OFS='\t' '{print $0, $2}' specific_targets_filt.txt > mosdepth_40_targets.bed

```
```
chr2_pat_hsa3	63820727	63820727
chr4_pat_hsa17x5	139673146	139673146
chr15_pat_hsa14	36212595	36212595
chr1_pat_hsa1	10233235	10233235
chr6_mat_hsa7	11390729	11390729
chr7_pat_hsa8	5143130	5143130
chr19_pat_hsa5x17	24137428	24137428
chr1_pat_hsa1	51992760	51992760
chr2_pat_hsa3	104430425	104430425
chr3_pat_hsa4	121313448	121313448
chr4_pat_hsa17x5	181243267	181243267
chr5_mat_hsa6	61734828	61734828
chr6_mat_hsa7	163795044	163795044
chr7_pat_hsa8	39191000	39191000
chr7_pat_hsa8	45067101	45067101
chr7_pat_hsa8	63184015	63184015
chr8_pat_hsa10	86733664	86733664
chr8_pat_hsa10	145958666	145958666
chr9_pat_hsa11	44720414	44720414
chr9_pat_hsa11	70173757	70173757
chr10_mat_hsa12	20613085	20613085
chr13_pat_hsa9	35940167	35940167
chr17_mat_hsa18	82026819	82026819
chr18_pat_hsa16	103784314	103784314
chr1_pat_hsa1	180923230	180923230
chr4_pat_hsa17x5	169976139	169976139
chr8_pat_hsa10	26547585	26547585
chr9_pat_hsa11  134720610	134720610
chr10_mat_hsa12	107788170	107788170
chr18_pat_hsa16	74778095	74778095
chr4_pat_hsa17x5	49217303	49217303
chr5_mat_hsa6	93647643	93647643
chr9_pat_hsa11	27699017	27699017
chr10_mat_hsa12	64331685	64331685
chr10_mat_hsa12	75066868	75066868
chr16_pat_hsa15	55694687	55694687
chr16_pat_hsa15	79107645	79107645
chr18_pat_hsa16	111163844	111163844
chr19_pat_hsa5x17	18421770	18421770
chr22_mat_hsa21	25719657	25719657
```

Now run `target_mosdepth_p1.sh`
```
#!/bin/bash
#SBATCH --job-name=mos
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=mos.%j.out # output file
#SBATCH --error=mos.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/bam_96_list_simplified.txt | sed -n ${ind}p)

mosdepth -n ${indiv_name} --by mosdepth_40_targets.bed /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/processed_bams/${indiv_name} && touch ${indiv_name}.mosdepth.done 
```
And submit `sbatch --array=1-96 target_mosdepth_p1.sh`
