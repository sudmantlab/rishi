Code for running [mpcrselect](https://github.com/ellieearmstrong/mPCRselect) pipeline
Make directory and prepare working environment
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/05_mpcrselect && cd /global/scratch/users/rdekayne/gorilla_census/05_mpcrselect

#https://github.com/ellieearmstrong/mPCRselect/tree/main
#install dependencies

#install nextflow
module load openjdk/17.0.8.1_1-gcc-11.4.0
curl -s https://get.nextflow.io | bash

#install pipeline
nextflow pull ellieearmstrong/mPCRselect -r main

#Clone the NGS-PrimerPlex repository
git clone https://github.com/aakechin/NGS-PrimerPlex

#Install the Python dependencies: 
cd NGS-PrimerPlex/
bash install_for_linux.sh 
chmod +x NGS_primerplex.py
#Modify the shebang line of the NGS_primerplex.py script: Change the first line from #!/usr/bin/python3 to #!/usr/bin/env python3
#Move the NGS_primerplex.py script to a directory in your PATH variable
cp NGS_primerplex.py /global/scratch/users/rdekayne/gorilla_census/

git clone https://github.com/ellieearmstrong/mPCRselect
cp mPCRselect/nextflow.config test.config

conda create -p /global/scratch/users/rdekayne/envs/mpcrselectenv
conda activate /global/scratch/users/rdekayne/envs/mpcrselectenv
conda install bioconda::bwa=0.7.17 
conda install conda-forge::libzlib=1.2.13 
conda install conda-forge::ruby=3.2.2 
conda install bioconda::vcftools=0.1.16 
conda install bioconda::plink2=2.00a5.10 
conda install bioconda::bcftools=1.18 
conda install conda-forge::gsl=2.7 
conda install conda-forge::gzip 
conda install conda-forge::gawk 
conda install conda-forge::r-tidyverse=1.3.1 
conda install conda-forge::r-caret 
conda install conda-forge::r-ggplot2=3.4.0
conda install bioconda::bedtools=2.31.0 
conda install -c conda-forge r-base=4.2.3
conda install bioconda::primer3-py
conda install conda-forge::biopython
conda install bioconda::pysam
conda install conda-forge::xlrd
conda install conda-forge::xlsxwriter
pip install networkx==1.11
#had to downgrade python:   python                pkgs/main::python-3.9.20-he870216_1 --> conda-forge::python-3.7.12-hf930737_100_cpython 
conda install python=3.7

git clone https://github.com/campanam/baitstools
cd baitstools
gem build baitstools.gemspec
gem install baitstools-1.8.1.gem

##check dependency versions match github
```

Produce sample csv file
```
Sample,Population
ind01_Maisha,mont1
ind02_Tuck,mont1
ind03_Turimaso,mont1
ind04_Umurimo,mont1
ind05_Imfura,mont1
ind06_Kaboko,mont1
ind07_Zirikana,mont1
ind17_Bwiruka,mont2
ind19_Katungi,mont2
ind20_Nyamunwa,mont2
ind21_Semehe,mont2
```

We will be editing one config file for the run - `test.config`

Copy the corilla genome and index here:
```
cp /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta ./gorilla.ref.fasta
cp /global/scratch/users/rdekayne/gorilla_census/data/genomes/mGorGor1.pri.cur.20231122.fasta.fai ./gorilla.ref.fasta.fai
```

We need to reheader the fasta and index file using this - `reheader_fasta.sh`
```
#!/bin/bash
#SBATCH --job-name=rehead
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=rehead.%j.out # output file
#SBATCH --error=rehead.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

sed -i 's/chr1_pat_hsa1/chr1pathsa1/g' ./gorilla.ref.fasta;
sed -i 's/chr2_pat_hsa3/chr2pathsa3/g' ./gorilla.ref.fasta;
sed -i 's/chr3_pat_hsa4/chr3pathsa4/g' ./gorilla.ref.fasta;
sed -i 's/chr4_pat_hsa17x5/chr4pathsa17x5/g' ./gorilla.ref.fasta;
sed -i 's/chr5_mat_hsa6/chr5mathsa6/g' ./gorilla.ref.fasta;
sed -i 's/chr6_mat_hsa7/chr6mathsa7/g' ./gorilla.ref.fasta;
sed -i 's/chr7_pat_hsa8/chr7pathsa8/g' ./gorilla.ref.fasta;
sed -i 's/chr8_pat_hsa10/chr8pathsa10/g' ./gorilla.ref.fasta;
sed -i 's/chr9_pat_hsa11/chr9pathsa11/g' ./gorilla.ref.fasta;
sed -i 's/chr10_mat_hsa12/chr10mathsa12/g' ./gorilla.ref.fasta;
sed -i 's/chr11_mat_hsa2b/chr11mathsa2b/g' ./gorilla.ref.fasta;
sed -i 's/chr12_pat_hsa2a/chr12pathsa2a/g' ./gorilla.ref.fasta;
sed -i 's/chr13_pat_hsa9/chr13pathsa9/g' ./gorilla.ref.fasta;
sed -i 's/chr14_pat_hsa13/chr14pathsa13/g' ./gorilla.ref.fasta;
sed -i 's/chr15_pat_hsa14/chr15pathsa14/g' ./gorilla.ref.fasta;
sed -i 's/chr16_pat_hsa15/chr16pathsa15/g' ./gorilla.ref.fasta;
sed -i 's/chr17_mat_hsa18/chr17mathsa18/g' ./gorilla.ref.fasta;
sed -i 's/chr18_pat_hsa16/chr18pathsa16/g' ./gorilla.ref.fasta;
sed -i 's/chr19_pat_hsa5x17/chr19pathsa5x17/g' ./gorilla.ref.fasta;
sed -i 's/chr20_mat_hsa19/chr20mathsa19/g' ./gorilla.ref.fasta;
sed -i 's/chr21_pat_hsa20/chr21pathsa20/g' ./gorilla.ref.fasta;
sed -i 's/chr22_mat_hsa21/chr22mathsa21/g' ./gorilla.ref.fasta;
sed -i 's/chr23_mat_hsa22/chr23mathsa22/g' ./gorilla.ref.fasta;
done

sed -i 's/chr1_pat_hsa1/chr1pathsa1/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr2_pat_hsa3/chr2pathsa3/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr3_pat_hsa4/chr3pathsa4/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr4_pat_hsa17x5/chr4pathsa17x5/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr5_mat_hsa6/chr5mathsa6/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr6_mat_hsa7/chr6mathsa7/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr7_pat_hsa8/chr7pathsa8/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr8_pat_hsa10/chr8pathsa10/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr9_pat_hsa11/chr9pathsa11/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr10_mat_hsa12/chr10mathsa12/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr11_mat_hsa2b/chr11mathsa2b/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr12_pat_hsa2a/chr12pathsa2a/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr13_pat_hsa9/chr13pathsa9/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr14_pat_hsa13/chr14pathsa13/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr15_pat_hsa14/chr15pathsa14/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr16_pat_hsa15/chr16pathsa15/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr17_mat_hsa18/chr17mathsa18/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr18_pat_hsa16/chr18pathsa16/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr19_pat_hsa5x17/chr19pathsa5x17/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr20_mat_hsa19/chr20mathsa19/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr21_pat_hsa20/chr21pathsa20/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr22_mat_hsa21/chr22mathsa21/g' ./gorilla.ref.fasta.fai;
sed -i 's/chr23_mat_hsa22/chr23mathsa22/g' ./gorilla.ref.fasta.fai;
done
```
And run - `sbatch reheader_fasta.sh`

Now prepare mpcrselect submission config and submission script to change paths etc.

Make an output directory for the output
```
mkdir run_one_output
```

And then run the pipeline - `submit_run_one_mpcrselect.sh `
```
#!/bin/bash
#SBATCH --job-name=mpcr
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=genomicdata_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=mpcr.%j.out
#SBATCH --error=mpcr.%j.err

module load openjdk/17.0.8.1_1-gcc-11.4.0

nextflow run ellieearmstrong/mPCRselect -r main -c run_one.config
```
Run - `sbatch submit_run_one_mpcrselect.sh`
