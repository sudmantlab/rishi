Here we will analyse the output of CAT which effectively expands the annotation from one of the genomes in our whole genome alignment (in our case the honeycomb rockfish) to the other taxa

As a result we get .gff3 files for each of the 55 haplotypes

Now we want to count the genes and look for overlaps of genes across the gene sets

Prepare the env
```
cd /global/scratch/users/rdekayne/rockfish55/assemblies
wget -i https://public.gi.ucsc.edu/~pnhebbar/rockfish_annotations/

mkdir -p /global/scratch/users/rdekayne/rockfish55/08_ANNOTATIONS && cd /global/scratch/users/rdekayne/rockfish55/08_ANNOTATIONS

conda create -p /global/scratch/users/rdekayne/envs/gff
conda activate /global/scratch/users/rdekayne/envs/gff
conda install bioconda::gffutils
```

Now run the following python script which will analyse these annotation files `analyse_gff3s.py`
```
import gffutils
import glob
import os

# Step 1: Find all .gff3 files in the current directory
gff_files = glob.glob("*.gff3")

# Step 2: Create a DB for each GFF and collect gene IDs
gene_sets = []
gene_counts = {}     # Store gene count per file
temp_db_paths = []   # Track temp DB files to delete later

for i, gff_file in enumerate(gff_files):
    db_path = f"temp_db_{i}.sqlite"
    temp_db_paths.append(db_path)

    # Create DB from GFF file
    gffutils.create_db(
        gff_file,
        dbfn=db_path,
        force=True,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True
    )

    # Open the database
    db = gffutils.FeatureDB(db_path)

    # Collect gene IDs into a set
    gene_ids = set()
    for gene in db.features_of_type('gene'):
        if 'ID' in gene.attributes:
            gene_ids.add(gene.attributes['ID'][0])

    gene_sets.append(gene_ids)
    gene_counts[gff_file] = len(gene_ids)

# Step 3: Find shared genes (intersection of all sets)
shared_genes = set.intersection(*gene_sets) if gene_sets else set()

# Step 4: Write shared genes to output file
with open("shared_genes.txt", "w") as f:
    for gene_id in sorted(shared_genes):
        f.write(gene_id + "\n")

# Step 5: Write gene counts per file
with open("gene_counts.txt", "w") as f:
    for gff_file in gff_files:
        f.write(f"{gff_file}: {gene_counts[gff_file]} genes\n")

# Step 6: Remove temporary SQLite DB files
for db_path in temp_db_paths:
    try:
        os.remove(db_path)
    except OSError as e:
        print(f"⚠️ Error deleting {db_path}: {e}")

# Final summary
print("\n📊 Gene counts per GFF3 file:")
for gff_file in gff_files:
    print(f"{gff_file}: {gene_counts[gff_file]} genes")

print(f"\nProcessed {len(gff_files)} GFF3 files.")
print(f"Found {len(shared_genes)} shared gene IDs.")
print("Output written to 'shared_genes.txt' and 'gene_counts.txt'")
```

Now we will run this from a shell script `run_analyse_gffs.sh`
```
#!/bin/bash
#SBATCH --job-name=gff3
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=gff3%j.out # output file
#SBATCH --error=gff3%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

python analyse_gff3s.py
touch FINISHED.txt
```

And run
```
sbatch run_analyse_gffs.sh
```
