Mapping of reads from Jan 2025 2-step PCR test

Make directories and load conda env
```
mkdir /global/scratch/users/rdekayne/gorilla_census/09_mapping && cd /global/scratch/users/rdekayne/gorilla_census/09_mapping
conda activate /global/scratch/users/rdekayne/envs/mapping
mkdir /global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams
```
Map data to reference genome - `09.1_mapping_map_p1.sh`
```
#!/bin/bash
#SBATCH --job-name=bwa_map
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=bwa_map.%j.out
#SBATCH --error=bwa_map.%j.err

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/sample_prefixes.txt | sed -n ${ind}p)

# call bwa
bwa-mem2 mem -t 24 /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/${indiv_name}_L001_R1_001.out.fastq.gz /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/${indiv_name}_L001_R2_001.out.fastq.gz | samtools sort -@24 -o /global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/${indiv_name}.raw.bam && touch ${indiv_name}.mapping.done
```
And run - `sbatch --array=1-39 09.1_mapping_map_p1.sh`

Tidy up output files
```

mkdir all_err_out
mv *.err all_err_out/
mv *.out all_err_out/

rm *.done
```
Load java for processing
```
module load java/22.0.1
```
Now do bam processing `09.2_mapping_process_p1.sh`
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
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/sample_prefixes.txt | sed -n ${ind}p)

#PART1
picard FixMateInformation I=/global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/${indiv_name}.raw.bam VALIDATION_STRINGENCY=LENIENT OUTPUT=/global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/tmp_bam/${indiv_name}.raw.bam

sambamba sort /global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/tmp_bam/${indiv_name}.raw.bam -o /global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/tmp_bam/${indiv_name}.sorted.raw.bam -t 24 -m 50GB

picard MarkDuplicates INPUT=/global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/tmp_bam/${indiv_name}.sorted.raw.bam OUTPUT=/global/scratch/users/rdekayne/gorilla_census/09_mapping/processed_bams/${indiv_name}.sorted.dup.bam METRICS_FILE=/global/scratch/users/rdekayne/gorilla_census/09_mapping/processed_bams/${indiv_name}.sorted.dup.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024

sambamba index /global/scratch/users/rdekayne/gorilla_census/09_mapping/processed_bams/${indiv_name}.sorted.dup.bam

#now remove all intermediate files
rm /global/scratch/users/rdekayne/gorilla_census/09_mapping/raw_bams/tmp_bam/${indiv_name}*

touch ${indiv_name}.processing.done
```
And submit - `sbatch --array=1-39%8 09.2_mapping_process_p1.sh`
Tidy output
```
mkdir all_err_out
mv *.err all_err_out/
mv *.out all_err_out/
```
Make a list of bam files which we will use later
```
ls processed_bams/*.bam > bam_39_list.txt
```
Now calculate mapping depth with mosdepth `09.3_mapping_mosdepth_p1.sh `
```
#!/bin/bash
#SBATCH --job-name=mos
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --output=mos.%j.out # output file
#SBATCH --error=mos.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat bam_39_list.txt | sed -n ${ind}p)

mosdepth -n ${indiv_name} /global/scratch/users/rdekayne/gorilla_census/09_mapping/processed_bams/${indiv_name} && touch ${indiv_name}.mosdepth.done 
samtools flagstat ${indiv_name} && touch ${indiv_name}.flagstst.done
```
Submit this `sbatch --array=1-39 09.3_mapping_mosdepth_p1.sh`
