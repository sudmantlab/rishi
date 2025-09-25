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

Load java - `module load java/22.0.1`
Then we are going to process the bam files including fixing mate information, marking duplicates, and sorting and indexing using a mixture of tools with `01.4_mapping_process_p1.sh `
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

#move to correct directory
cd /global/scratch/users/rdekayne/gorilla_census/01_mapping

#set up array process using samples
ind=${SLURM_ARRAY_TASK_ID}
#cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt
cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/sample_list2.txt | sed -n ${ind}p > ${SLURM_ARRAY_TASK_ID}.sample_list.txt

VAR1=$(cut -d ' ' -f1 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)
VAR2=$(cut -d ' ' -f2 ${SLURM_ARRAY_TASK_ID}.sample_list.txt)

#PART1
picard FixMateInformation I=/global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/${VAR1}_${VAR2}.raw.bam VALIDATION_STRINGENCY=LENIENT OUTPUT=/global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/tmp_bam/${VAR1}_${VAR2}.raw.bam

sambamba sort/global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/tmp_bam/${VAR1}_${VAR2}.raw.bam -o /global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/tmp_bam/${VAR1}_${VAR2}.sorted.raw.bam -t 24 -m 50GB

picard MarkDuplicates INPUT=/global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/tmp_bam/${VAR1}_${VAR2}.sorted.raw.bam OUTPUT=/global/scratch/users/rdekayne/gorilla_census/01_mapping/processed_bams/${VAR1}_${VAR2}.sorted.dup.bam METRICS_FILE=/global/scratch/users/rdekayne/gorilla_census/01_mapping/processed_bams/${VAR1}_${VAR2}.sorted.dup.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024

sambamba index /global/scratch/users/rdekayne/gorilla_census/01_mapping/processed_bams/${VAR1}_${VAR2}.sorted.dup.bam 

#now remove all intermediate files
rm /global/scratch/users/rdekayne/gorilla_census/01_mapping/raw_bams/tmp_bam/${VAR1}_${VAR2}*
rm ${SLURM_ARRAY_TASK_ID}.sample_list.txt

touch ${VAR1}_${VAR2}.processing.done
```
Submit that `sbatch --array=1-131%8 01.4_mapping_process_p1.sh`

Now move directory
`mkdir -p /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams`
  
We are now going to merge bam files for individuals with multiple bams - `01.5_mapping_merge_p1.sh`
```
#!/bin/bash
#SBATCH --job-name=bwa_merg
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3_bigmem
#SBATCH --output=bwa_merg.%j.out # output file
#SBATCH --error=bwa_merg.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/merge.list | sed -n ${ind}p)

samtools merge -o /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/${indiv_name}_1file.bam /global/scratch/users/rdekayne/gorilla_census/01_mapping/processed_bams/${indiv_name}*.bam

sambamba index /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/${indiv_name}_1file.bam

touch ${indiv_name}.mege.done
```
Submit `sbatch --array=1-9 01.5_mapping_merge_p1.sh`

Rename the files so they are more intuitive (i.e. by individual) and then make a list of them to use for genotyping
```
cp ind01_Maisha_ERR668423.sorted.dup.bam ../indiv_bams/ind01_Maisha_1file.bam
cp ind01_Maisha_ERR668423.sorted.dup.bam.bai ../indiv_bams/ind01_Maisha_1file.bam.bai
cp ind01_Maisha_ERR668423.sorted.dup.txt ../indiv_bams/ind01_Maisha_1file.txt
cp ind03_Turimaso_ERR668425.sorted.dup.bam ../indiv_bams/ind03_Turimaso_1file.bam
cp ind03_Turimaso_ERR668425.sorted.dup.bam.bai ../indiv_bams/ind03_Turimaso_1file.bam.bai
cp ind03_Turimaso_ERR668425.sorted.dup.txt ../indiv_bams/ind03_Turimaso_1file.txt
cp ind04_Umurimo_ERR668424.sorted.dup.bam ../indiv_bams/ind04_Umurimo_1file.bam
cp ind04_Umurimo_ERR668424.sorted.dup.bam.bai ../indiv_bams/ind04_Umurimo_1file.bam.bai
cp ind04_Umurimo_ERR668424.sorted.dup.txt ../indiv_bams/ind04_Umurimo_1file.txt
cp ind08_Dunia_ERR668428.sorted.dup.bam ../indiv_bams/ind08_Dunia_1file.bam
cp ind08_Dunia_ERR668428.sorted.dup.bam.bai ../indiv_bams/ind08_Dunia_1file.bam.bai
cp ind08_Dunia_ERR668428.sorted.dup.txt ../indiv_bams/ind08_Dunia_1file.txt
cp ind10_Pinga_ERR668427.sorted.dup.bam ../indiv_bams/ind10_Pinga_1file.bam
cp ind10_Pinga_ERR668427.sorted.dup.bam.bai ../indiv_bams/ind10_Pinga_1file.bam.bai
cp ind10_Pinga_ERR668427.sorted.dup.txt ../indiv_bams/ind10_Pinga_1file.txt
cp ind11_Serufuli_ERR668429.sorted.dup.bam ../indiv_bams/ind11_Serufuli_1file.bam
cp ind11_Serufuli_ERR668429.sorted.dup.bam.bai ../indiv_bams/ind11_Serufuli_1file.bam.bai
cp ind11_Serufuli_ERR668429.sorted.dup.txt ../indiv_bams/ind11_Serufuli_1file.txt
cp ind12_Tumani_ERR668426.sorted.dup.bam ../indiv_bams/ind12_Tumani_1file.bam
cp ind12_Tumani_ERR668426.sorted.dup.bam.bai ../indiv_bams/ind12_Tumani_1file.bam.bai
cp ind12_Tumani_ERR668426.sorted.dup.txt ../indiv_bams/ind12_Tumani_1file.txt

cp ind17_Bwiruka_ERR2300765.sorted.dup.bam ../indiv_bams/ind17_Bwiruka_1file.bam
cp ind17_Bwiruka_ERR2300765.sorted.dup.bam.bai ../indiv_bams/ind17_Bwiruka_1file.bam.bai
cp ind17_Bwiruka_ERR2300765.sorted.dup.txt ../indiv_bams/ind17_Bwiruka_1file.txt
cp ind19_Katungi_ERR2300763.sorted.dup.bam ../indiv_bams/ind19_Katungi_1file.bam
cp ind19_Katungi_ERR2300763.sorted.dup.bam.bai ../indiv_bams/ind19_Katungi_1file.bam.bai
cp ind19_Katungi_ERR2300763.sorted.dup.txt ../indiv_bams/ind19_Katungi_1file.txt
cp ind20_Nyamunwa_ERR2300764.sorted.dup.bam ../indiv_bams/ind20_Nyamunwa_1file.bam
cp ind20_Nyamunwa_ERR2300764.sorted.dup.bam.bai ../indiv_bams/ind20_Nyamunwa_1file.bam.bai
cp ind20_Nyamunwa_ERR2300764.sorted.dup.txt ../indiv_bams/ind20_Nyamunwa_1file.txt
cp ind21_Semehe_ERR2300766.sorted.dup.bam ../indiv_bams/ind21_Semehe_1file.bam
cp ind21_Semehe_ERR2300766.sorted.dup.bam.bai ../indiv_bams/ind21_Semehe_1file.bam.bai
cp ind21_Semehe_ERR2300766.sorted.dup.txt ../indiv_bams/ind21_Semehe_1file.txt

#make list for genotyping
cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams && ls *.bam > bam_1file.list
```

Calculate depth across bam files - `01.6_mapping_mosdepth_p1.sh `
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
indiv_name=$(cat bam_1file.list | sed -n ${ind}p)

mosdepth -n ${indiv_name} /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/${indiv_name} && touch ${indiv_name}.mosdepth.done 
samtools flagstat ${indiv_name} && touch ${indiv_name}.flagstst.done
```
Submit `sbatch --array=1-16 01.6_mapping_mosdepth_p1.sh`

Make a separate list for just runnning the mountain gorillas `bam_2file.list`
```
ind17_Bwiruka_1file.bam
ind19_Katungi_1file.bam
ind20_Nyamunwa_1file.bam
ind21_Semehe_1file.bam`
```

And run mosdepth on just these indviduals too `01.6_mapping_mosdepth_p2.sh`
```
#!/bin/bash
#SBATCH --job-name=mos
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=mos.%j.out # output file
#SBATCH --error=mos.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams

ind=${SLURM_ARRAY_TASK_ID}
indiv_name=$(cat bam_2file.list | sed -n ${ind}p)

mosdepth -n ${indiv_name} /global/scratch/users/rdekayne/gorilla_census/01_mapping/indiv_bams/${indiv_name} && touch ${indiv_name}.mosdepth.done 
samtools flagstat ${indiv_name} && touch ${indiv_name}.flagstst.done
```
Submit this too `sbatch --array=1-4 01.6_mapping_mosdepth_p2.sh`







