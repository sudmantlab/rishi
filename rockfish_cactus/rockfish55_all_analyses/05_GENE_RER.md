We are now going to run RER converge on genes

To do this we will pull out exon alignments from our .hal file - if these exon alignments have all 55 taxa we will then cat them together by gene and make a tree from each gene alignment

First we need to prepare the environment 
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/05_GENE_RER && cd /global/scratch/users/rdekayne/rockfish55/05_GENE_RER
cd /global/scratch/users/rdekayne/rockfish/03_CACTUS/cactus
source cactus_env/bin/activate
mkdir -p /global/scratch/users/rdekayne/rockfish55/05_GENE_RER/test_tmp
```

Now we will process the annotation file for honeycomb rockfish 

First take honeycomb gff (honeycomb.gff) and filter out gene and exon lines
```
grep -v "#" honeycomb.gff > honeycomb_nohash.gff
awk '$3 == "gene" || $3 == "exon"' honeycomb_nohash.gff > gene_exon_honey.gff
awk '$3 == "gene" { match($0, /gene=([^;]+)/, arr); if (arr[1]) print arr[1] }' gene_exon_honey.gff | sort -u > sorted_unique_gene_names.txt
```

We again need this to be in batches of genes < 999 in length to run with sarray so split this file
```
split -l 999 -d sorted_unique_gene_names.txt cut_sorted_genes
mkdir -p gene_exon_maf_files
```

Then use this template scriptn `run_hal2maf_genes00.sh` to run each batch of genes one by one
```
#!/bin/bash
#SBATCH --job-name=cactus_100kb
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_gene%j.out # output file
#SBATCH --error=cactus_gene%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24 

export TMPDIR=/global/scratch/users/rdekayne/rockfish55/05_GENE_RER/test_tmp
#get number
numb=${SLURM_ARRAY_TASK_ID}

#extract gene name
## CHANGE THIS LINE FOR EACH BATCH OF 999 GENES
gene=$(cat cut_sorted_genes00 | sed -n ${numb}p)

#get exons from gene_exon_honey.gff
grep ";gene=${gene};" gene_exon_honey.gff > tmp_${gene}_exons.txt

# Remove the 'gene' line and keep only exons
grep -v -P "\tgene\t" tmp_${gene}_exons.txt > tmp_${gene}_exon_only.txt

# Count the number of exons
num_exons=$(wc -l < tmp_${gene}_exon_only.txt)

# Loop through each exon
for i in $(seq 1 $num_exons); do
    # Zero-pad exon number (e.g. 001, 002, ...)
    exon_num=$(printf "%03d" $i)

    # Extract the exon line
    exon_line=$(sed -n "${i}p" tmp_${gene}_exon_only.txt)

    # Get the chromosome, start, and end (columns 1, 4, 5)
    chr=$(echo "$exon_line" | cut -f1)
    start=$(echo "$exon_line" | cut -f4)
    end=$(echo "$exon_line" | cut -f5)

    # Save to a new tmp file
    echo -e "${chr}\t${start}\t${end}" > tmp_${gene}_exon${exon_num}.bed

    singularity exec ../../rockfish/03_CACTUS/cactus_v2.9.3.sif hal2maf ../03_CACTUS/55_rockfish.hal ./gene_exon_maf_files/${gene}_exon${exon_num}.maf --refGenome honeycomb-vgp.1 --noAncestors --noDupes --unique --refTargets tmp_${gene}_exon${exon_num}.bed

done

#tidy output
rm tmp_${gene}_exon*
```
And submit this file
```
sbatch --array=1%20 run_hal2maf_genes00.sh
sbatch --array=1-999%100 run_hal2maf_genes00.sh
```

Then we want to convert these maf files to a msa file for each exon, make sure all taxa are present in the msa and if so combine into a single gene msa
```
mkdir -p exon_msas
mkdir -p gene_msas
cp ../04_HAL_TO_TREES/55_sed_subs.sed . 

conda activate /global/scratch/users/rdekayne/envs/maf_manipulation
export TMPDIR=/global/scratch/users/rdekayne/rockfish55/05_GENE_RER/test_tmp
```

Here we have a big script `convert_maf2fasta00.sh` that does all of this
```
#!/bin/bash
#SBATCH --job-name=cactus_test
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=cactus_test%j.out # output file
#SBATCH --error=cactus_test%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

export TMPDIR=/global/scratch/users/rdekayne/rockfish55/05_GENE_RER/test_tmp

#get number from job
numb=${SLURM_ARRAY_TASK_ID}

#extract gene name
## CHANGE THIS LINE FOR EACH BATCH OF 999 GENES
gene=$(cat cut_sorted_genes00 | sed -n ${numb}p)

#get exons from gene_exon_honey.gff
grep ";gene=${gene};" gene_exon_honey.gff > tmp_${gene}_exons.txt

#remove the 'gene' line and keep only exons
grep -v -P "\tgene\t" tmp_${gene}_exons.txt > tmp_${gene}_exon_only.txt

#count number of exons
num_exons=$(wc -l < tmp_${gene}_exon_only.txt)

#then make a list of .maf files for that gene
ls ./gene_exon_maf_files/${gene}_exon*.maf > full_${gene}_exon_list.txt

#count how many there are (some exons may have not been convered to mafs e.g. no data)
num_exon_files=$(wc -l < full_${gene}_exon_list.txt)

#loop through each exon
for i in $(seq 1 $num_exon_files); do
    #select the correct exon
    exon_maf_filename=$(cat full_${gene}_exon_list.txt | sed -n ${i}p)

    #make a bunch of sed substitutions for naming changes
    sed -i -f 55_sed_subs.sed "${exon_maf_filename}"

    #now use msa view for conversions - we only want to keep those specific sequences
    msa_view ${exon_maf_filename} --seqs korean,SebastesborealisSEB-8_shortraker_1,SebastesborealisSEB-8_shortraker_2,SebastesmelanostictusSEB-4_blackspotted_1,SebastesmelanostictusSEB-4_blackspotted_2,7164-PS-0003_rougheye_1,7164-PS-0003_rougheye_2,SebastespolyspinisSEB-5_northern_1,SebastespolyspinisSEB-5_northern_2,7164-PS-0007_lightdusky_1,7164-PS-0007_lightdusky_2,SebastesalutusPOP6_popa_1,SebastesalutusPOP6_popa_2,7164-PS-0006_popb_1,7164-PS-0006_popb_2,acadian-vgp_1,acadian-vgp_2,acadian-cbp_1,acadian-cbp_2,SebastesserranoidesSEB-259_olive_1,SebastesserranoidesSEB-259_olive_2,7164-PS-0002_black_1,7164-PS-0002_black_2,7164-PS-0005_yellowtail_1,7164-PS-0005_yellowtail_2,7164-PS-0004_widow_1,7164-PS-0004_widow_2,widow-ccgp_1,widow-ccgp_2,SebastesmystinusSEB-254_blue_1,SebastesmystinusSEB-254_blue_2,m84185_240412_222911_s4_deacon_1,m84185_240412_222911_s4_deacon_2,SebastescrocotulusSEB-252_sunset_1,SebastescrocotulusSEB-252_sunset_2,SebastesminiatusSEB-256_vermilion_1,SebastesminiatusSEB-256_vermilion_2,bocaccio-ccgp_1,bocaccio-ccgp_2,Sebastesruberrimus12-Yellowye1a_1,Sebastesruberrimus12-Yellowye1a_2,honeycomb-vgp_1,honeycomb-vgp_2,SebastescarnatusSEB-258_gopher_1,SebastescarnatusSEB-258_gopher_2,7164-PS-0001_copper_1,7164-PS-0001_copper_2,7164-PS-0008_quillback_1,7164-PS-0008_quillback_2,SebastesalascanusSEB-2_shortspinethornyhead_1,SebastesalascanusSEB-2_shortspinethornyhead_2,Sebastesaleutianus_rougheye_kolora,Sebastespinniger_canary_kolora,Sebastesminiatus_vermilion_kolora,Sebastesrosaceus_rosy_kolora, -f -G 1 --unmask | seqkit seq -w 0 > ${exon_maf_filename}.fasta

    #tidy up the output specifically missing bases
    sed -e 's/\*/N/g' ${exon_maf_filename}.fasta > ${exon_maf_filename}.gap-replaced.fasta
    sed -i 's/> />/g' ${exon_maf_filename}.gap-replaced.fasta
    rm ${exon_maf_filename}.fasta
    mv ./gene_exon_maf_files/${gene}_exon*.gap-replaced.fasta ./exon_msas/

    #tidy output
    #rm tmp_${gene}_exon*
done

##in this part we are going to look at the gene level and concatenate info across exons
#get output file name
output="${gene}_full.fasta"

#make tmp directory for storing intermediate files
tmp_dir=$(mktemp -d)

#loop through all exon fasta files matching the gene prefix set above
for fasta in ./exon_msas/${gene}_exon*.fasta; do

    #count number of FASTA headers (i.e., number of sequences) - we need this to be 50 to have full data
    num_headers=$(grep -c "^>" "$fasta")

    #conditional to include only include files with exactly 55 sequences
    if [ "$num_headers" -eq 55 ]; then
        #generate md5sum file of headers to check pre-concat that the order is the same
        grep ">" "$fasta" > "${fasta}_md5_tmp.txt"

        #another conditional to check the md5 with a known order
        if [ "$(md5sum ${fasta}_md5_tmp.txt | awk '{print $1}')" = "$(md5sum md5_fasta_headers.txt | awk '{print $1}')" ]; then
            echo "Files have the same MD5 checksum."
            echo "Including $fasta (55 sequences)"
            awk -v dir="$tmp_dir" '
            /^>/ {header=$0; gsub(/^>/, "", header); next}
            {
                print >> dir "/" header
            }' "$fasta"
        else
            echo "Files are different."
        fi
    else
        echo "Skipping $fasta (has $num_headers sequences)"
    fi
    #rm ${fasta}_md5_tmp.txt
done

#concatenate the sequences by taxon into the final output
> "$output"  # Create or clear the output file
for taxon_file in "$tmp_dir"/*; do
    taxon=$(basename "$taxon_file")
    echo ">$taxon" >> "$output"
    tr -d '\n' < "$taxon_file" >> "$output"
    echo >> "$output"
done

#now clean up the output
mv "$output" ./gene_msas/"$output"
rm -r "$tmp_dir"
rm full_${gene}_exon_list.txt 
echo "Combined alignment written to $output"
```

And run this as a huge array in batches of 999 genes
``` 
sbatch --array=1-10%20 convert_maf2fasta00.sh
sbatch --array=1-999%100 convert_maf2fasta00.sh
```
create directory for tiny files (small)
```
cd /global/scratch/users/rdekayne/rockfish55/05_GENE_RER/gene_msas
mkdir -p small_files
```

Loop through each fasta file and copy small files to `small_files`
```
for file in *.fasta; do
    # Check if it's a regular file
    if [ -f "$file" ]; then
        # Get the file size in bytes (Linux)
        size=$(stat -c%s "$file")

        # If size is less than or equal to 5, copy to small_files/
        if [ "$size" -le 20000 ]; then
            mv "$file" small_files/
        fi
    fi
done
```

ADD VARIABLE SITES ETC.
