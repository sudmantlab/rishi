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
Then prepare input bed file for mosdepth `mosdepth_40_targets.bed`
```
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
indiv_name=$(cat bam_96_list.txt | sed -n ${ind}p)

mosdepth -n ${indiv_name} --by mosdepth_40_targets.bed /path/to/bams/${indiv_name} && touch ${indiv_name}.mosdepth.done 
```
##run 
sbatch --array=1-16 01.6_mapping_mosdepth_p1.sh 
