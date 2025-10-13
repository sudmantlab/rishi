In this section we will calculate [buscos](https://busco.ezlab.org/) across the assemblies

First make our directories and load a caonda environment to work in
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/01_BUSCOS && cd /global/scratch/users/rdekayne/rockfish55/01_BUSCOS
mkdir -p /global/scratch/users/rdekayne/rockfish55/01_BUSCOS/output_busco_phylogenomics

#start by calculating buscos across the assemblies starting with the 19 that nicolas sent:
conda create -p /global/scratch/users/rdekayne/envs/busco
conda activate /global/scratch/users/rdekayne/envs/busco
conda install -c conda-forge -c bioconda busco=5.8.2
```

Now run the busco step with `55_rockfish_buscos.sh`
```
#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=busco%j.out # output file
#SBATCH --error=busco%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

numb=${SLURM_ARRAY_TASK_ID}
file_name=$(cat /global/scratch/users/rdekayne/rockfish55/assemblies/rockfish55_assembly.list | sed -n ${numb}p)

#now busco command
busco -c 24 -i /global/scratch/users/rdekayne/rockfish55/assemblies/${file_name} -o ${file_name} --out_path output_busco_phylogenomics/ -l s/global/scratch/users/rdekayne/rockfish55/01_BUSCOS/busco_downloads/lineages/actinopterygii_odb10/ -m genome

touch ./${file_name}.done
```

And then run with `sbatch --array=1-55%10 55_rockfish_buscos.sh`

Now parse the output
```
touch filenames.txt
touch scores.txt
for file in *.out
do
grep "Input file" ${file} >> filenames.txt
grep "C:" ${file} >> scores.txt
done

paste filenames.txt scores.txt > full_busco_output_summary.txt
```

