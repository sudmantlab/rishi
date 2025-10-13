Here we will extract information from our .hal cactus alignment

In particular we want to pull out 100kb windows along the length of the genome

First prepare the environment and make folders
```
#make folders
mkdir -p /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES

#load env
cd /global/scratch/users/rdekayne/rockfish/03_CACTUS/cactus
source cactus_env/bin/activate
cd /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES
mkdir -p mkdir whole_hal_to_maf && cd whole_hal_to_maf
mkdir test_tmp
```

Now convert the entire .hal to .maf files with `conv_hal_to_maf_honey.sh` oriented/polarised on the honeycomb rockfish (for which we have an annotation

```
#!/bin/bash
#SBATCH --job-name=cactus_hal_to_maf
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_test%j.out # output file
#SBATCH --error=cactus_test%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

singularity exec ../../../rockfish/03_CACTUS/cactus_v2.9.3.sif cactus-hal2maf ./rockfish55_js ../../03_CACTUS/55_rockfish.hal 55_rockfish_honey.maf --refGenome honeycomb-vgp.1 --noAncestors --dupeMode single --filterGapCausingDupes --chunkSize 1000000 --batchMemory 25G --maxMemory 50G
```
And submit with `sbatch conv_hal_to_maf_honey.sh`

Now to get our windowed analysis we need to first make our windows

We use samtools to index the genome and bedtools to make windows from our file
```
cd /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES
cp ../assemblies/honeycomb-vgp_1.fasta .
conda activate /global/scratch/users/rdekayne/envs/trees
samtools faidx honeycomb-vgp_1.fasta
cut -f1,2 honeycomb-vgp_1.fasta.fai > sizes.genome
bedtools makewindows -g sizes.genome -w 100000 > 100kb.txt
```
Now check for one chromosome
```
#for chr1
grep "NC_051269.1" 100kb.txt | wc -l
##441
```

And the total number of window
```
wc -l 100kb.txt
##8084 100kb.txt
```

To run on our cluster we cant have an array >999 jobs so we need to split by chromosome
```
#split the 100kb text file by chromosome:
awk '{print >> ($1 "_chr.txt")}' 100kb.txt
ls NC_* > all_NC_chr.txt
mkdir small_chrs
mv NW_* small_chrs/
```
Now load env and make directories
```
cd /global/scratch/users/rdekayne/rockfish/03_CACTUS/cactus
source cactus_env/bin/activate
mkdir -p /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/maf_run100 && cd /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/maf_run100
mkdir bed_files
mkdir maf_files
```

Now HAL --> MAF with `run_hal2maf_windows.sh`
```

#!/bin/bash
#SBATCH --job-name=cactus_100kb
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_test%j.out # output file
#SBATCH --error=cactus_test%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

export TMPDIR=/global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/test_tmp

# Get task ID from SLURM array
numb=${SLURM_ARRAY_TASK_ID}

# Base directory
BASE_DIR="/global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES"

# File containing list of window files
FILE_LIST="${BASE_DIR}/all_NC_chr.txt"

# Get the filename for this task
chr_filename=$(sed -n "${numb}p" "$FILE_LIST")
chr_file="${BASE_DIR}/${chr_filename}"

# Check that the file exists
if [[ ! -f "$chr_file" ]]; then
    echo "Error: File $chr_file does not exist."
    exit 1
fi

# Read each line (window) from the selected chromosome file
line_number=0
while IFS= read -r window; do
    ((line_number++))

    # Create a .bed file for this line
    bed_file="chr_${numb}_l_${line_number}.bed"
    echo "$window" > "$bed_file"

    # Output MAF file
    maf_file="chr_${numb}_l_${line_number}.maf"

    # Run hal2maf
    singularity exec ../../../rockfish/03_CACTUS/cactus_v2.9.3.sif \
        hal2maf ../../03_CACTUS/55_rockfish.hal "$maf_file" \
        --refGenome honeycomb-vgp.1 \
        --noAncestors \
        --noDupes \
        --unique \
        --refTargets "$bed_file"

done < "$chr_file"

mv *.bed ./bed_files/
mv *.maf ./maf_files/
```
And run this
```
sbatch --array=1-1%4 run_hal2maf_windows.sh
sbatch --array=2-24%1 run_hal2maf_windows.sh
```

Now MAF --> MSA

Prepare the environment
```
conda activate /global/scratch/users/rdekayne/envs/maf_manipulation
mkdir all_MSAs
mkdir 55_MSAs
mkdir fasta_files
```

Since we have a sed step to replace a lot of names prepare a .sed file `55_sed_subs.sed` to save compute
```
s/SebastesborealisSEB-8_shortraker.1/SebastesborealisSEB-8_shortraker_1/g
s/SebastesborealisSEB-8_shortraker.2/SebastesborealisSEB-8_shortraker_2/g
s/SebastesmelanostictusSEB-4_blackspotted.1/SebastesmelanostictusSEB-4_blackspotted_1/g
s/SebastesmelanostictusSEB-4_blackspotted.2/SebastesmelanostictusSEB-4_blackspotted_2/g
s/7164-PS-0003_rougheye.1/7164-PS-0003_rougheye_1/g
s/7164-PS-0003_rougheye.2/7164-PS-0003_rougheye_2/g
s/SebastespolyspinisSEB-5_northern.1/SebastespolyspinisSEB-5_northern_1/g
s/SebastespolyspinisSEB-5_northern.2/SebastespolyspinisSEB-5_northern_2/g
s/7164-PS-0007_lightdusky.1/7164-PS-0007_lightdusky_1/g
s/7164-PS-0007_lightdusky.2/7164-PS-0007_lightdusky_2/g
s/SebastesalutusPOP6_popa.1/SebastesalutusPOP6_popa_1/g
s/SebastesalutusPOP6_popa.2/SebastesalutusPOP6_popa_2/g
s/7164-PS-0006_popb.1/7164-PS-0006_popb_1/g
s/7164-PS-0006_popb.2/7164-PS-0006_popb_2/g
s/acadian-vgp.1/acadian-vgp_1/g
s/acadian-vgp.2/acadian-vgp_2/g
s/acadian-cbp.1/acadian-cbp_1/g
s/acadian-cbp.2/acadian-cbp_2/g
s/SebastesserranoidesSEB-259_olive.1/SebastesserranoidesSEB-259_olive_1/g
s/SebastesserranoidesSEB-259_olive.2/SebastesserranoidesSEB-259_olive_2/g
s/7164-PS-0002_black.1/7164-PS-0002_black_1/g
s/7164-PS-0002_black.2/7164-PS-0002_black_2/g
s/7164-PS-0005_yellowtail.1/7164-PS-0005_yellowtail_1/g
s/7164-PS-0005_yellowtail.2/7164-PS-0005_yellowtail_2/g
s/7164-PS-0004_widow.1/7164-PS-0004_widow_1/g
s/7164-PS-0004_widow.2/7164-PS-0004_widow_2/g
s/widow-ccgp.1/widow-ccgp_1/g
s/widow-ccgp.2/widow-ccgp_2/g
s/SebastesmystinusSEB-254_blue.1/SebastesmystinusSEB-254_blue_1/g
s/SebastesmystinusSEB-254_blue.2/SebastesmystinusSEB-254_blue_2/g
s/m84185_240412_222911_s4_deacon.1/m84185_240412_222911_s4_deacon_1/g
s/m84185_240412_222911_s4_deacon.2/m84185_240412_222911_s4_deacon_2/g
s/SebastescrocotulusSEB-252_sunset.1/SebastescrocotulusSEB-252_sunset_1/g
s/SebastescrocotulusSEB-252_sunset.2/SebastescrocotulusSEB-252_sunset_2/g
s/SebastesminiatusSEB-256_vermilion.1/SebastesminiatusSEB-256_vermilion_1/g
s/SebastesminiatusSEB-256_vermilion.2/SebastesminiatusSEB-256_vermilion_2/g
s/bocaccio-ccgp.1/bocaccio-ccgp_1/g
s/bocaccio-ccgp.2/bocaccio-ccgp_2/g
s/Sebastesruberrimus12-Yellowye1a.1/Sebastesruberrimus12-Yellowye1a_1/g
s/Sebastesruberrimus12-Yellowye1a.2/Sebastesruberrimus12-Yellowye1a_2/g
s/honeycomb-vgp.1/honeycomb-vgp_1/g
s/honeycomb-vgp.2/honeycomb-vgp_2/g
s/SebastescarnatusSEB-258_gopher.1/SebastescarnatusSEB-258_gopher_1/g
s/SebastescarnatusSEB-258_gopher.2/SebastescarnatusSEB-258_gopher_2/g
s/7164-PS-0001_copper.1/7164-PS-0001_copper_1/g
s/7164-PS-0001_copper.2/7164-PS-0001_copper_2/g
s/7164-PS-0008_quillback.1/7164-PS-0008_quillback_1/g
s/7164-PS-0008_quillback.2/7164-PS-0008_quillback_2/g
s/SebastesalascanusSEB-2_shortspinethornyhead.1/ebastesalascanusSEB-2_shortspinethornyhead_1/g
s/SebastesalascanusSEB-2_shortspinethornyhead.2/ebastesalascanusSEB-2_shortspinethornyhead_2/g
```

Now run the conversion script `convert_maf2fasta.sh`
```
#!/bin/bash
#SBATCH --job-name=cactus_mafmsa
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_mafmsa%j.out # output file
#SBATCH --error=cactus_mafmsa%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

export TMPDIR=/global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/test_tmp

# Get task ID from SLURM array
numb=${SLURM_ARRAY_TASK_ID}

# Base directory
BASE_DIR="/global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES"

# File containing list of window files
FILE_LIST="${BASE_DIR}/all_NC_chr.txt"

# Get the filename for this task
chr_filename=$(sed -n "${numb}p" "$FILE_LIST")
chr_file="${BASE_DIR}/${chr_filename}"

# Check that the file exists
if [[ ! -f "$chr_file" ]]; then
    echo "Error: File $chr_file does not exist."
    exit 1
fi

# Read each line (window) from the selected chromosome file
line_number=0
while IFS= read -r window; do
    ((line_number++))

    maf_file="chr_${numb}_l_${line_number}.maf"
    fasta_file="chr_${numb}_l_${line_number}.fasta"
    gap_replaced_fasta="chr_${numb}_l_${line_number}.gap-replaced.fasta"

    #sed subs
    sed -i -f 55_sed_subs.sed ./maf_files/"$maf_file"

    #run msa_view
    msa_view ./maf_files/"$maf_file" --seqs korean,SebastesborealisSEB-8_shortraker_1,SebastesborealisSEB-8_shortraker_2,SebastesmelanostictusSEB-4_blackspotted_1,SebastesmelanostictusSEB-4_blackspotted_2,7164-PS-0003_rougheye_1,7164-PS-0003_rougheye_2,SebastespolyspinisSEB-5_northern_1,SebastespolyspinisSEB-5_northern_2,7164-PS-0007_lightdusky_1,7164-PS-0007_lightdusky_2,SebastesalutusPOP6_popa_1,SebastesalutusPOP6_popa_2,7164-PS-0006_popb_1,7164-PS-0006_popb_2,acadian-vgp_1,acadian-vgp_2,acadian-cbp_1,acadian-cbp_2,SebastesserranoidesSEB-259_olive_1,SebastesserranoidesSEB-259_olive_2,7164-PS-0002_black_1,7164-PS-0002_black_2,7164-PS-0005_yellowtail_1,7164-PS-0005_yellowtail_2,7164-PS-0004_widow_1,7164-PS-0004_widow_2,widow-ccgp_1,widow-ccgp_2,SebastesmystinusSEB-254_blue_1,SebastesmystinusSEB-254_blue_2,m84185_240412_222911_s4_deacon_1,m84185_240412_222911_s4_deacon_2,SebastescrocotulusSEB-252_sunset_1,SebastescrocotulusSEB-252_sunset_2,SebastesminiatusSEB-256_vermilion_1,SebastesminiatusSEB-256_vermilion_2,bocaccio-ccgp_1,bocaccio-ccgp_2,Sebastesruberrimus12-Yellowye1a_1,Sebastesruberrimus12-Yellowye1a_2,honeycomb-vgp_1,honeycomb-vgp_2,SebastescarnatusSEB-258_gopher_1,SebastescarnatusSEB-258_gopher_2,7164-PS-0001_copper_1,7164-PS-0001_copper_2,7164-PS-0008_quillback_1,7164-PS-0008_quillback_2,SebastesalascanusSEB-2_shortspinethornyhead_1,SebastesalascanusSEB-2_shortspinethornyhead_2,Sebastesaleutianus_rougheye_kolora,Sebastespinniger_canary_kolora,Sebastesminiatus_vermilion_kolora,Sebastesrosaceus_rosy_kolora, -f -G 1 --unmask | seqkit seq -w 0 > ./fasta_files/"$fasta_file"

    sed -e 's/\*/N/g' ./fasta_files/"$fasta_file" > "$gap_replaced_fasta"
    sed -i 's/> />/g' "$gap_replaced_fasta"
    mv "$gap_replaced_fasta" ./all_MSAs/"$gap_replaced_fasta"

done < "$chr_file"
```

And run this with 
```
sbatch --array=1-1%4 convert_maf2fasta.sh
sbatch --array=2-24%1 convert_maf2fasta.sh
```

Now we want to filter the msa outputs

First make a directory
```
mkdir 50_taxa_msas
```

Now run `check_msas.sh` to make sure there are 55 headers in there i.e. no missing taxa and if no missing copy to our `55_MSAs` folder
```
#!/bin/bash

# Output file
output_file="all_taxa_msas.txt"

# Clear the output file if it exists
> "$output_file"

# Loop through all .fasta files in the current directory
for file in ./all_MSAs/*.fasta; do
    # Count lines starting with '>'
    header_count=$(grep -c '^>' "$file")

    # If count is exactly 50, append filename to output
    if [ "$header_count" -eq 55 ]; then
        echo "$file" >> "$output_file"
        cp "$file" ./55_MSAs
    fi
done
```

Now to actually get phylogenetic trees from these windows I'll use pargenes

First load the conda env
```
conda activate /global/scratch/users/rdekayne/envs/trees
```

Then run pargenes on our no missing data folder with `pargenes.sh`
```
#!/bin/bash
#SBATCH --job-name=pargenes
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=pargenes%j.out # output file
#SBATCH --error=pargenes%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=56 # 2 CPUs per job

/global/scratch/users/rdekayne/rockfish/ParGenes/pargenes/pargenes.py -a /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/maf_run100/msa_files_100k/55_MSAs -o ./pargenes_output --cores 56 -d nt -m --job-failure-fatal --use-astral
```
And submit `sbatch pargenes.sh`

#now make a directory for the pargenes output since we will export these trees locally  and ptu the trees there
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/whole_genome_tree_windows && cd /global/scratch/users/rdekayne/rockfish55/04_HAL_TO_TREES/whole_genome_tree_windows
cp ../maf_run100/pargenes_output/mlsearch_run/results/*/*.bestTree .
```
