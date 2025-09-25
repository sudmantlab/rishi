De-Kayne: Code to download data to make gorilla SNP panel

Make directory
```
cd /global/scratch/users/rdekayne/
mkdir -p /global/scratch/users/rdekayne/envs/
mkdir -p /global/scratch/users/rdekayne/gorilla_census && cd /global/scratch/users/rdekayne/gorilla_census
mkdir -p /global/scratch/users/rdekayne/gorilla_census/data && cd /global/scratch/users/rdekayne/gorilla_census/data
```

Prepare conda environment for download
```
module load anaconda3
conda create -p /global/scratch/users/rdekayne/envs/sra
source activate base
conda activate /global/scratch/users/rdekayne/envs/sra
conda install bioconda::sra-tools
```

Download Gorilla genome we will use
```
#Genomes
mkdir -p /global/scratch/users/rdekayne/gorilla_census/data/genomes && cd /global/scratch/users/rdekayne/gorilla_census/data/genomes
#https://github.com/marbl/Primates?tab=readme-ov-file
#https://genomeark.s3.amazonaws.com/index.html?prefix=species/Gorilla_gorilla/mGorGor1/assembly_curated/
```

We are going to download files per individual/accession from SRA

`01_download_gorilla_data_01.sh`
```
#!/bin/bash
#SBATCH --job-name=sra
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl.out # output file
#SBATCH --error=dl.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

mkdir -p ind14_Mkubwa && cd ind14_Mkubwa
fasterq-dump --split-files SRX243452
touch SRX243452.done
fasterq-dump --split-files SRX243453
touch SRX243453.done
```

`01_download_gorilla_data_02.sh`
```
#!/bin/bash
#SBATCH --job-name=sra2
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl2.out # output file
#SBATCH --error=dl2.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

mkdir -p ind15_Kaisi && cd ind15_Kaisi
fasterq-dump --split-files SRX242685
touch SRX242685.done
fasterq-dump --split-files SRX242686
touch SRX242686.done
fasterq-dump --split-files SRX242687
touch SRX242687.done
fasterq-dump --split-files SRX242688
touch SRX242688.done
```

`01_download_gorilla_data_03.sh`
```
#!/bin/bash
#SBATCH --job-name=sra
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl3.out # output file
#SBATCH --error=dl3.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

mkdir -p ind16_Victoria && cd ind16_Victoria
fasterq-dump --split-files SRX243528
touch SRX243528.done
fasterq-dump --split-files SRX243529
touch SRX243529.done
fasterq-dump --split-files SRX243530
touch SRX243530.done
fasterq-dump --split-files SRX243531
touch SRX243531.done
fasterq-dump --split-files SRX243532
touch SRX243532.done
fasterq-dump --split-files SRX243533
touch SRX243533.done
```

`01_download_gorilla_data_04.sh`
```
#!/bin/bash
#SBATCH --job-name=wget04
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl4.out # output file
#SBATCH --error=dl4.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind02_Tuck && cd ind02_Tuck
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223795/ERR223795_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223735/ERR223735_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223759/ERR223759_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223741/ERR223741_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223711/ERR223711_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223771/ERR223771_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223717/ERR223717_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223729/ERR223729_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223801/ERR223801_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223807/ERR223807_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223747/ERR223747_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223753/ERR223753_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223729/ERR223729_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223777/ERR223777_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223765/ERR223765_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223807/ERR223807_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225697/ERR225697_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223783/ERR223783_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223771/ERR223771_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223801/ERR223801_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223789/ERR223789_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223717/ERR223717_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223747/ERR223747_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223765/ERR223765_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223723/ERR223723_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223783/ERR223783_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223723/ERR223723_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223759/ERR223759_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223741/ERR223741_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223735/ERR223735_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225697/ERR225697_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223789/ERR223789_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223711/ERR223711_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223753/ERR223753_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223777/ERR223777_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223795/ERR223795_2.fastq.gz
touch ERS168204_others.done
```

`01_download_gorilla_data_05.sh`
```
#!/bin/bash
#SBATCH --job-name=wget05
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl5.out # output file
#SBATCH --error=dl5.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind03_Turimaso && cd ind03_Turimaso
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668425/ERR668425_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668425/ERR668425_2.fastq.gz
touch ERR668425.done
```

`01_download_gorilla_data_06.sh`
```
#!/bin/bash
#SBATCH --job-name=wget06
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl6.out # output file
#SBATCH --error=dl6.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind05_Imfura && cd ind05_Imfura
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223776/ERR223776_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223806/ERR223806_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223746/ERR223746_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225696/ERR225696_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223740/ERR223740_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223782/ERR223782_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223710/ERR223710_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223734/ERR223734_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225696/ERR225696_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223746/ERR223746_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223752/ERR223752_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223800/ERR223800_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223764/ERR223764_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223770/ERR223770_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223776/ERR223776_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223722/ERR223722_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223794/ERR223794_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223728/ERR223728_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223758/ERR223758_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223716/ERR223716_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223740/ERR223740_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223788/ERR223788_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223716/ERR223716_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223758/ERR223758_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223764/ERR223764_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223722/ERR223722_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223794/ERR223794_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223728/ERR223728_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223770/ERR223770_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223710/ERR223710_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223800/ERR223800_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223752/ERR223752_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223806/ERR223806_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223788/ERR223788_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223734/ERR223734_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223782/ERR223782_1.fastq.gz
touch ERS168207_others.done
```

`01_download_gorilla_data_07.sh`
```
#!/bin/bash
#SBATCH --job-name=wet07
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl7.out # output file
#SBATCH --error=dl7.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind06_Kaboko && cd ind06_Kaboko
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223768/ERR223768_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223792/ERR223792_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223726/ERR223726_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223720/ERR223720_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223750/ERR223750_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223786/ERR223786_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223738/ERR223738_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223762/ERR223762_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223744/ERR223744_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225700/ERR225700_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223732/ERR223732_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223714/ERR223714_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223780/ERR223780_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223774/ERR223774_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223798/ERR223798_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223810/ERR223810_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223804/ERR223804_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223774/ERR223774_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223714/ERR223714_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223756/ERR223756_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223726/ERR223726_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223738/ERR223738_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223720/ERR223720_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223732/ERR223732_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223768/ERR223768_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223756/ERR223756_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223750/ERR223750_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223786/ERR223786_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223810/ERR223810_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223792/ERR223792_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223804/ERR223804_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223798/ERR223798_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223744/ERR223744_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223762/ERR223762_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223780/ERR223780_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225700/ERR225700_2.fastq.gz
touch ERS168410_others.done
```

`01_download_gorilla_data_08.sh`
```
#!/bin/bash
#SBATCH --job-name=wget08
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl8.out # output file
#SBATCH --error=dl8.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind07_Zirikana && cd ind07_Zirikana
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225695/ERR225695_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223787/ERR223787_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223745/ERR223745_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223727/ERR223727_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223793/ERR223793_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223733/ERR223733_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223751/ERR223751_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223775/ERR223775_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223769/ERR223769_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223805/ERR223805_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223751/ERR223751_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223745/ERR223745_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223793/ERR223793_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223799/ERR223799_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223709/ERR223709_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223757/ERR223757_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223739/ERR223739_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223715/ERR223715_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223763/ERR223763_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223733/ERR223733_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223739/ERR223739_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223721/ERR223721_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223805/ERR223805_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223775/ERR223775_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223781/ERR223781_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223799/ERR223799_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223757/ERR223757_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225695/ERR225695_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223781/ERR223781_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223787/ERR223787_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223727/ERR223727_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223769/ERR223769_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223715/ERR223715_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223763/ERR223763_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223721/ERR223721_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223709/ERR223709_1.fastq.gz
touch ERS168174_others.done
```

`01_download_gorilla_data_09.sh`
```
#!/bin/bash
#SBATCH --job-name=wget09
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl9.out # output file
#SBATCH --error=dl9.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind08_Dunia && cd ind08_Dunia
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668428/ERR668428_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668428/ERR668428_2.fastq.gz
touch ERS525621_others.done
```

`01_download_gorilla_data_10.sh`
```
#!/bin/bash
#SBATCH --job-name=wget10
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl10.out # output file
#SBATCH --error=dl10.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind09_Itebero && cd ind09_Itebero
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223784/ERR223784_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225698/ERR225698_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223742/ERR223742_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223754/ERR223754_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223718/ERR223718_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223808/ERR223808_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223778/ERR223778_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223796/ERR223796_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223736/ERR223736_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223760/ERR223760_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223790/ERR223790_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223778/ERR223778_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223766/ERR223766_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223808/ERR223808_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223748/ERR223748_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223784/ERR223784_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223802/ERR223802_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223712/ERR223712_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223760/ERR223760_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223742/ERR223742_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223724/ERR223724_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223790/ERR223790_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223730/ERR223730_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223766/ERR223766_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223754/ERR223754_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223748/ERR223748_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223796/ERR223796_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223772/ERR223772_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223802/ERR223802_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223724/ERR223724_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223730/ERR223730_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225698/ERR225698_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223772/ERR223772_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223718/ERR223718_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223736/ERR223736_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223712/ERR223712_1.fastq.gz
touch ERS168205_others.done
```

`01_download_gorilla_data_11.sh`
```
#!/bin/bash
#SBATCH --job-name=wget11
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl11.out # output file
#SBATCH --error=dl11.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind10_Pinga && cd ind10_Pinga
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668427/ERR668427_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668427/ERR668427_1.fastq.gz
touch ERS525620_others.done
```

`01_download_gorilla_data_12.sh`
```
#!/bin/bash
#SBATCH --job-name=wget12
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl12.out # output file
#SBATCH --error=dl12.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind13_Ntabwoba && cd ind13_Ntabwoba
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223809/ERR223809_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223779/ERR223779_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223761/ERR223761_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223737/ERR223737_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223791/ERR223791_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225699/ERR225699_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223713/ERR223713_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223743/ERR223743_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223773/ERR223773_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223779/ERR223779_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223803/ERR223803_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223749/ERR223749_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223767/ERR223767_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223719/ERR223719_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223713/ERR223713_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223731/ERR223731_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223737/ERR223737_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223785/ERR223785_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223755/ERR223755_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223725/ERR223725_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223797/ERR223797_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225699/ERR225699_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223719/ERR223719_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223803/ERR223803_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223761/ERR223761_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223731/ERR223731_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223785/ERR223785_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223773/ERR223773_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223755/ERR223755_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223767/ERR223767_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223809/ERR223809_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223743/ERR223743_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223749/ERR223749_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223725/ERR223725_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223797/ERR223797_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223791/ERR223791_1.fastq.gz
touch ERS168206_others.done
```

`01_download_gorilla_data_13.sh`
```
#!/bin/bash
#SBATCH --job-name=wget1
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl13.out # output file
#SBATCH --error=dl13.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind01_Maisha && cd ind01_Maisha

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668423/ERR668423_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668423/ERR668423_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668423/ERR668423.fastq.gz
touch ind01.done
```

`01_download_gorilla_data_14.sh`
```
#!/bin/bash
#SBATCH --job-name=wget2
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl14.out # output file
#SBATCH --error=dl14.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind04_Umurimo && cd ind04_Umurimo

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668424/ERR668424_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668424/ERR668424.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668424/ERR668424_2.fastq.gz
touch ind04.done
```

`01_download_gorilla_data_15.sh`
```
#!/bin/bash
#SBATCH --job-name=wget3
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl15.out # output file
#SBATCH --error=dl15.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind11_Serufuli && cd ind11_Serufuli
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668429/ERR668429_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668429/ERR668429.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668429/ERR668429_1.fastq.gz
touch ind11.done
```

`01_download_gorilla_data_16.sh`
```
#!/bin/bash
#SBATCH --job-name=wget4
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=dl16.out # output file
#SBATCH --error=dl16.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

mkdir -p ind12_Tumani && cd ind12_Tumani
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668426/ERR668426_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668426/ERR668426.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR668/ERR668426/ERR668426_1.fastq.gz
touch ind12.done
```

Now some extra scripts for downloading tricky files
```
#!/bin/bash
#SBATCH --job-name=BONUS1
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=BONUS1.out # output file
#SBATCH --error=BONUS1.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

cd ind02_Tuck
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223807/ERR223807_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223747/ERR223747_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223753/ERR223753_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223729/ERR223729_2.fastq.gz

cd ../ind05_Imfura
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223752/ERR223752_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223800/ERR223800_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223764/ERR223764_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223770/ERR223770_1.fastq.gz
```

```
#!/bin/bash
#SBATCH --job-name=BONUS2
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=BONUS2.out # output file
#SBATCH --error=BONUS2.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

cd ind06_Kaboko
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223762/ERR223762_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223744/ERR223744_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225700/ERR225700_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223732/ERR223732_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223714/ERR223714_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223780/ERR223780_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223774/ERR223774_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223798/ERR223798_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223810/ERR223810_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223804/ERR223804_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223774/ERR223774_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223714/ERR223714_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223756/ERR223756_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223726/ERR223726_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223738/ERR223738_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223720/ERR223720_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223732/ERR223732_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223768/ERR223768_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223756/ERR223756_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223750/ERR223750_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223786/ERR223786_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223810/ERR223810_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223792/ERR223792_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223804/ERR223804_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223798/ERR223798_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223744/ERR223744_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223762/ERR223762_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223780/ERR223780_2.fastq.gz
```


```
#!/bin/bash
#SBATCH --job-name=BONUS3
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=BONUS3.out # output file
#SBATCH --error=BONUS3.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

cd ind07_Zirikana
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223775/ERR223775_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223769/ERR223769_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223805/ERR223805_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223751/ERR223751_1.fastq.gz

cd ../ind13_Ntabwoba
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223743/ERR223743_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223773/ERR223773_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223779/ERR223779_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223803/ERR223803_1.fastq.gz
```

```
#!/bin/bash
#SBATCH --job-name=BONUS4
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --output=BONUS4.out # output file
#SBATCH --error=BONUS4.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit of 4GB

cd ind09_Itebero
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223796/ERR223796_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223736/ERR223736_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223760/ERR223760_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223790/ERR223790_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223778/ERR223778_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223766/ERR223766_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223808/ERR223808_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223748/ERR223748_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223784/ERR223784_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223802/ERR223802_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223712/ERR223712_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223760/ERR223760_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223742/ERR223742_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223724/ERR223724_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223790/ERR223790_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223730/ERR223730_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223766/ERR223766_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223754/ERR223754_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223748/ERR223748_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223796/ERR223796_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/ERR223772/ERR223772_2.fastq.gz
```

And delete files that are not paired
```
##deleted unpaired files:
rm ind01_Maisha/ERR668423.fastq.gz
rm ind04_Umurimo/ERR668424.fastq.gz
rm ind12_Tumani/ERR668426.fastq.gz
rm ind15_Kaisi/SRR747657.fastq.gz
rm ind15_Kaisi/SRR747658.fastq.gz
```

And download additional mountain gorilla accessions



