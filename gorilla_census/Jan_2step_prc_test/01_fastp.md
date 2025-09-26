Code to do QC of reads from Jan 2025 test of 2-step pcr primers for amplifying 156 SNP loci

Set up env and make directories
```
mkdir /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered
cd /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered

conda activate /global/scratch/users/rdekayne/envs/fastp
```

Prepare sample lists from raw data:
```
ls ../Fastq/*.fastq.gz > full_sample.txt
#get one sample ID per individual so remove R2 files
grep "R1" full_sample.txt > full_sampleR1.txt
#remove file extention
sed -i 's/_L001_R1_001.fastq.gz//g' full_sampleR1.txt
sed -i 's/..\/Fastq\///g' full_sampleR1.txt
grep -v "PCR1null_S37_L001" full_sampleR1.txt > sample_prefixes.txt
#check how many
wc -l sample_prefixes.txt
```
Results in `39 sample_prefixes.txt`

Read QC with fastp
```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=fastq1.%j.out # output file
#SBATCH --error=fastq1.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=16 # 2 CPUs per job

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/sample_prefixes.txt | sed -n ${ind}p)

fastp -i /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq/${indiv_name}_L001_R1_001.fastq.gz -I /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq/${indiv_name}_L001_R2_001.fastq.gz --thread 16 -o /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/${indiv_name}_L001_R1_001.out.fastq.gz -O /global/scratch/users/rdekayne/gorilla_census/Jan_2025_library_test_UCB/Fastq_fastp_filtered/${indiv_name}_L001_R2_001.out.fastq.gz -j ${indiv_name}.json -h ${indiv_name}.html && touch ${indiv_name}.done
```
And run with `sbatch --array=1-39 08.1_mapping_qc_p1.sh `

Make folders for output
```
mkdir all_jsons
mkdir all_html

mv *.json all_jsons/
mv *.html all_html/

mkdir all_err_out
mv *.err all_err_out/
mv *.out all_err_out/
```

And run multiqc to get a single report
```
eval "$(/global/scratch/users/scott_ferguson/modules/miniconda3/bin/conda shell.zsh hook)"
mamba activate multiqc
cd all_jsons/
multiqc --fn_as_s_name . 
```
