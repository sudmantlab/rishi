Code for mapping of downloaded SRA gorilla data  

Make directories and prepare conda env(s)
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/01_mapping && cd /global/scratch/users/rdekayne/gorilla_census/01_mapping
mkdir /global/scratch/users/rdekayne/gorilla_census/01_mapping/fastp_fastqs

conda create -p /global/scratch/users/rdekayne/envs/fastp
conda activate /global/scratch/users/rdekayne/envs/fastp
conda install bioconda::fastp
```

Do QC of reads using fastp - `01.1_mapping_qc_p1.sh `
```
#!/bin/bash
#SBATCH --job-name=fastp_test
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --qos=savio_lowprio
#SBATCH --output=fastq1.%j.out # output file
#SBATCH --error=fastq1.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=16 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping

ind=${SLURM_ARRAY_TASK_ID}
#cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt
cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list2.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt

VAR1=$(cut -d ' ' -f1 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)
VAR2=$(cut -d ' ' -f2 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)

echo $VAR1
echo $VAR2

fastp -i /global/scratch/users/rdekayne/gorilla_census/data/${VAR1}/${VAR2}_1.fastq.gz -I /global/scratch/users/rdekayne/gorilla_census/data/${VAR1}/${VAR2}_2.fastq.gz --thread 16 -o /global/scratch/users/rdekayne/gorilla_census/01_mapping/fastp_fastqs/${VAR1}_${VAR2}_1.out.fastq.gz -O /global/scratch/users/rdekayne/gorilla_census/01_mapping/fastp_fastqs/${VAR1}_${VAR2}_2.out.fastq.gz -j ${VAR1}_${VAR2}.json -h ${VAR1}_${VAR2}.html && touch ${VAR1}_${VAR2}.done

rm ${SLURM_ARRAY_TASK_ID}.sample_list.txt
```
and run using `sbatch --array=1-131 01.1_mapping_qc_p1.sh`
  
Combine into one report with multiqc
``` 
eval "$(/global/scratch/users/scott_ferguson/modules/miniconda3/bin/conda shell.zsh hook)"
mamba activate multiqc

mkdir -p /global/scratch/users/rdekayne/gorilla_census/01_mapping/html_output && cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/html_output
mv *.html /global/scratch/users/rdekayne/gorilla_census/01_mapping/html_output

mkdir -p /global/scratch/users/rdekayne/gorilla_census/01_mapping/json_output && cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/json_output
mv *.json /global/scratch/users/rdekayne/gorilla_census/01_mapping/json_output

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/json_output
multiqc --fn_as_s_name .
```

Conda env:
```
conda create -p /global/scratch/users/rdekayne/envs/mapping
conda activate /global/scratch/users/rdekayne/envs/mapping
conda install bioconda::samtools
conda install bioconda::sambamba
conda install bioconda::bwa-mem2
conda install bioconda::bcftools
conda install bioconda::picard
conda install bioconda::mosdepth
```

Index the gorilla reference - `01.2_mapping_index_p1.sh `
```
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --qos=savio_lowprio
#SBATCH --output=index.out # output file
#SBATCH --error=index.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer

cd /global/scratch/users/rdekayne/gorilla_census/data/genomes
bwa-mem2 index mGorGor1.pri.cur.20231122.fasta
```

```
mkdir /global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams
```

Now map the reads to the indexed genome - `01.3_mapping_map_p1.sh`
```
#!/bin/bash
#SBATCH --job-name=bwa_map
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=bwa_map.%j.out
#SBATCH --error=bwa_map.%j.err
#SBATCH -w n0067.savio4

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping

ind=${SLURM_ARRAY_TASK_ID}
#cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt
cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list2.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt

VAR1=$(cut -d ' ' -f1 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)
VAR2=$(cut -d ' ' -f2 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)

# call bwa
bwa-mem2 mem -t 24 /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta /global/scratch/users/rdekayne/gorilla_census/01_mapping/fastp_fastqs/${VAR1}_${VAR2}_1.out.fastq.gz /global/scratch/users/rdekayne/gorilla_census/01_mapping/fastp_fastqs/${VAR1}_${VAR2}_2.out.fastq.gz | samtools sort -@24 -o /global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/${VAR1}_${VAR2}.raw.bam && touch ${VAR1}_${VAR2}.mapping.done
```
and run with `sbatch --array=1-131 01.3_mapping_map_p1.sh `








