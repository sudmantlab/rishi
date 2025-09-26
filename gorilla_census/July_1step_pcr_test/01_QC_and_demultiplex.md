Code outlining QC of reads for gorilla 1step PCR test  

Load conda env and make directory
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/ && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/

conda activate /global/scratch/users/rdekayne/envs/fastp
```
Now run QC `QC_fastp.sh`
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

fastp -i 18032D-186-01_S234_L007_R1_001.fastq.gz -I 18032D-186-01_S234_L007_R2_001.fastq.gz --thread 16 -o 18032D-186-01_S234_L007_R1_001.out.fastq.gz -O 18032D-186-01_S234_L007_R2_001.out.fastq.gz -j OUTPUT.json -h OUTPUT.html && touch fastp.done
```
and run `sbatch QC_fastp.sh`

Now going to demultiplex the reads. Make directory
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/S1_folder && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/S1_folder
```
Now we need to demultiplex our reads - this is a bit complicated because illumina barcodes have ligated in both orientations  

This requires a directory `/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCBbarcode_primer_files/` which contains text files named F_01-F_12 and R_01-R_08 containing lists of primer sequence e.g. `F_01` starts with:
```
AGACTATGAGGAACTGAAACTTACCAGATTAC
AGACTATGACGCGGTTGGTGTAGTAT
AGACTATGGCTCCAAGCTGTGCATTT
```
A copy of this directory with the barcodes can be found [HERE](https://github.com/sudmantlab/rishi/tree/main/gorilla_census/July_1step_pcr_test/barcode_primer_files)

Now run the R script `/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/make_primer_list.R` to prepare forwards and reverse primer lists

```
#De-Kayne 2025
#produce list of barcodes for every combination of F and R primers for 40 loci

#set directory
setwd("/global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/")


fwd_lst <- c("F_01", "F_02", "F_03", "F_04", "F_05", "F_06", "F_07", "F_08", "F_09", "F_10", "F_11", "F_12")
rev_lst <- c("R_01", "R_02", "R_03", "R_04", "R_05", "R_06", "R_07", "R_08")

full_list_df <- as.data.frame(c())
full_list_df$name <- c()
full_list_df$primers <- c()

for(i in 1:length(fwd_lst)){
  for(j in 1:length(rev_lst)){
    tmp_DF <- as.data.frame(c(1:40))
    fileF <- read.csv(paste("./barcode_primer_files/", fwd_lst[i], sep = ""), head = F)
    fileR <- read.csv(paste("./barcode_primer_files/", rev_lst[j], sep = ""), head = F)
    tmp_DF$primer_F <- fileF$V1
    tmp_DF$primer_R <- fileR$V1
    tmp_DF$dash <- "-"
    tmp_DF$PRIMER <- paste(tmp_DF$primer_F, tmp_DF$dash, tmp_DF$primer_R, sep = "")
    tmp_DF$NAME <- paste(fwd_lst[i], "_", rev_lst[j], "_", tmp_DF$`c(1:40)`, sep = "")
    exp_DF <- as.data.frame(cbind(tmp_DF$NAME, tmp_DF$PRIMER))
    colnames(exp_DF) <- c("name", "primers")     
    full_list_df <- as.data.frame(rbind(full_list_df, exp_DF))
  }
}

write.table(full_list_df, file = "all_primers_for_demultiplexing.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
full_list_df_FWD <- full_list_df

#now just flip order of fwd and rev
full_list_df <- as.data.frame(c())
full_list_df$name <- c()
full_list_df$primers <- c()

for(i in 1:length(fwd_lst)){
  for(j in 1:length(rev_lst)){
    tmp_DF <- as.data.frame(c(1:40))
    fileF <- read.csv(paste("./barcode_primer_files/", fwd_lst[i], sep = ""), head = F)
    fileR <- read.csv(paste("./barcode_primer_files/", rev_lst[j], sep = ""), head = F)
    tmp_DF$primer_F <- fileF$V1
    tmp_DF$primer_R <- fileR$V1
    tmp_DF$dash <- "-"
    tmp_DF$PRIMER <- paste(tmp_DF$primer_R, tmp_DF$dash, tmp_DF$primer_F, sep = "")
    tmp_DF$NAME <- paste(fwd_lst[i], "_", rev_lst[j], "_", tmp_DF$`c(1:40)`, sep = "")
    exp_DF <- as.data.frame(cbind(tmp_DF$NAME, tmp_DF$PRIMER))
    colnames(exp_DF) <- c("name", "primers")     
    full_list_df <- as.data.frame(rbind(full_list_df, exp_DF))
  }
}

write.table(full_list_df, file = "all_primers_for_demultiplexing_flip.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
Now split both the forwards file and reverse file in 4 so we can run it effectively (there are many combinations)
```
split -n 4 all_primers_for_demultiplexing.txt all_primers_for_demultiplexing_
split -n 4 all_primers_for_demultiplexing_flip.txt all_primers_for_demultiplexing_REVCOMP_
```
Make a new directory for them and copy them here for demultiplexing
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex && cd /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex
cp ../all_primers_for_demultiplexing_* .
cp ../all_primers_for_demultiplexing_REVCOMP_* .
```
Now make new directories we will demultiplex into
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/demux
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/demux_rev
```
Demultiplex the forwards reads `demultiplex_all_fwd.sh`
```
#!/bin/bash
#SBATCH --job-name=demult
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_gene%j.out # output file
#SBATCH --error=cactus_gene%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_aa ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux/%_R1.fastq.gz ./demux/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_ab ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux/%_R1.fastq.gz ./demux/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_ac ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux/%_R1.fastq.gz ./demux/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_ad ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux/%_R1.fastq.gz ./demux/%_R2.fastq.gz -m 1
```
and run `sbatch demultiplex_all_fwd.sh`  

Now demultiplex the reverse reads `demultiplex_all_revcomp.sh`
```
#!/bin/bash
#SBATCH --job-name=demult
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_gene%j.out # output file
#SBATCH --error=cactus_gene%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_REVCOMP_aa ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux_rev/%_R1.fastq.gz ./demux_rev/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_REVCOMP_ab ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux_rev/%_R1.fastq.gz ./demux_rev/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_REVCOMP_ac ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux_rev/%_R1.fastq.gz ./demux_rev/%_R2.fastq.gz -m 1

../fastq-multx/fastq-multx -B all_primers_for_demultiplexing_REVCOMP_ad ../18032D-186-01_S234_L007_R1_001.fastq.gz ../18032D-186-01_S234_L007_R2_001.fastq.gz -o ./demux_rev/%_R1.fastq.gz ./demux_rev/%_R2.fastq.gz -m 1
```
and submit `sbatch demultiplex_all_revcomp.sh`  

This produces many files for each individual barcode set and locus
```
ls demux/* | wc -l
```
:7682
```
ls demux_rev/* | wc -l
```
:7682
  
Now we will concatenate all the reads per primer set i.e. per individual/replicate  

We need forwards and reverse orientation together  

Make directory  
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/July_2025_library_test_UCB/all_demultiplex/combined_indiv_fastas
```
And then make and run `concatentate_reads.sh`
```
#!/bin/bash

# Enable nullglob so that unmatched globs expand to nothing instead of themselves
shopt -s nullglob

# Create output directory if it doesn't exist
mkdir -p ./combined_indiv_fastas

# Loop over F_01 to F_12
for i in $(printf "%02d\n" {1..12}); do
  F="F_$i"

  # Loop over R_01 to R_08
  for j in $(printf "%02d\n" {1..8}); do
    R="R_$j"
    prefix="${F}_${R}"

    echo "Processing ${prefix}..."

    # Define output filenames
    out_R1="./combined_indiv_fastas/concat_${prefix}_R1.fastq.gz"
    out_R2="./combined_indiv_fastas/concat_${prefix}_R2.fastq.gz"

    # Collect matching files
    R1_files=(demux/${prefix}_*R1.fastq.gz demux_rev/${prefix}_*R1.fastq.gz)
    R2_files=(demux/${prefix}_*R2.fastq.gz demux_rev/${prefix}_*R2.fastq.gz)

    # Check if files exist before concatenating
    if [ ${#R1_files[@]} -eq 0 ]; then
      echo "⚠️  No R1 files found for ${prefix}, skipping..."
    else
      cat "${R1_files[@]}" > "$out_R1"
      echo "✅ Created $out_R1"
    fi

    if [ ${#R2_files[@]} -eq 0 ]; then
      echo "⚠️  No R2 files found for ${prefix}, skipping..."
    else
      cat "${R2_files[@]}" > "$out_R2"
      echo "✅ Created $out_R2"
    fi

    echo ""
  done
done
```
Make the file executible and execute
```
chmod +x concatentate_reads.sh 
./concatentate_reads.sh
```
Now check our output:
```
ls combined_indiv_fastas/* | wc -l
```
192 which = 2x96 which is our full plate of primers as we should have
