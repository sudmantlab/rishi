Make some new directories for storing the assemblies we will be analysing

```
mkdir -p /global/scratch/users/rdekayne/rockfish55 && cd /global/scratch/users/rdekayne/rockfish55
mkdir -p /global/scratch/users/rdekayne/rockfish55/assemblies && cd /global/scratch/users/rdekayne/rockfish55/assemblies
```

Now copy the assemblies to our new directory using `copy_assemblies.sh`
```
#copy_assemblies.sh
#!/bin/bash
#SBATCH --job-name=copy_assemblies
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=copy_assemblies%j.out # output file
#SBATCH --error=copy_assemblies%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=4 # 2 CPUs per job

#copy assemblies
#hap resolved
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0001.fa.masked ./7164-PS-0001_copper_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0001.fa.masked ./7164-PS-0001_copper_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0002.fa.masked ./7164-PS-0002_black_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0002.fa.masked ./7164-PS-0002_black_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0003.fa.masked ./7164-PS-0003_rougheye_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0003.fa.masked ./7164-PS-0003_rougheye_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0004.fa.masked ./7164-PS-0004_widow_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0004.fa.masked ./7164-PS-0004_widow_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0005.fa.masked ./7164-PS-0005_yellowtail_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0005.fa.masked ./7164-PS-0005_yellowtail_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0006.fa.masked ./7164-PS-0006_popb_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0006.fa.masked ./7164-PS-0006_popb_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0007.fa.masked ./7164-PS-0007_lightdusky_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0007.fa.masked ./7164-PS-0007_lightdusky_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/7164-PS-0008.fa.masked ./7164-PS-0008_quillback_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/7164-PS-0008.fa.masked ./7164-PS-0008_quillback_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/m84185_240412_222911_s4.fa.masked ./m84185_240412_222911_s4_deacon_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/m84185_240412_222911_s4.fa.masked ./m84185_240412_222911_s4_deacon_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesalascanusSEB-2.fa.masked ./SebastesalascanusSEB-2_shortspinethornyhead_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesalascanusSEB-2.fa.masked ./SebastesalascanusSEB-2_shortspinethornyhead_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesalutusPOP6.fa.masked ./SebastesalutusPOP6_popa_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesalutusPOP6.fa.masked ./SebastesalutusPOP6_popa_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesborealisSEB-8.fa.masked ./SebastesborealisSEB-8_shortraker_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesborealisSEB-8.fa.masked ./SebastesborealisSEB-8_shortraker_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastescarnatusSEB-258.fa.masked ./SebastescarnatusSEB-258_gopher_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastescarnatusSEB-258.fa.masked ./SebastescarnatusSEB-258_gopher_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastescrocotulusSEB-252.fa.masked ./SebastescrocotulusSEB-252_sunset_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastescrocotulusSEB-252.fa.masked ./SebastescrocotulusSEB-252_sunset_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesmelanostictusSEB-4.fa.masked ./SebastesmelanostictusSEB-4_blackspotted_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesmelanostictusSEB-4.fa.masked ./SebastesmelanostictusSEB-4_blackspotted_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesminiatusSEB-256.fa.masked ./SebastesminiatusSEB-256_vermilion_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesminiatusSEB-256.fa.masked ./SebastesminiatusSEB-256_vermilion_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesmystinusSEB-254.fa.masked ./SebastesmystinusSEB-254_blue_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesmystinusSEB-254.fa.masked ./SebastesmystinusSEB-254_blue_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastespolyspinisSEB-5.fa.masked ./SebastespolyspinisSEB-5_northern_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastespolyspinisSEB-5.fa.masked ./SebastespolyspinisSEB-5_northern_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/Sebastesruberrimus12-Yellowye1a.fa.masked ./Sebastesruberrimus12-Yellowye1a_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/Sebastesruberrimus12-Yellowye1a.fa.masked ./Sebastesruberrimus12-Yellowye1a_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap1/SebastesserranoidesSEB-259.fa.masked ./SebastesserranoidesSEB-259_olive_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/SebastesserranoidesSEB-259.fa.masked ./SebastesserranoidesSEB-259_olive_2.fasta
#primary/alt
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/acadian-cbp.fa.masked ./acadian-cbp_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/acadian-cbp.fa.masked ./acadian-cbp_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/acadian-vgp.fa.masked ./acadian-vgp_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/acadian-vgp.fa.masked ./acadian-vgp_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/bocaccio-ccgp.fa.masked ./bocaccio-ccgp_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/bocaccio-ccgp.fa.masked ./bocaccio-ccgp_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/honeycomb-vgp.fa.masked ./honeycomb-vgp_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/honeycomb-vgp.fa.masked ./honeycomb-vgp_2.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/widow-ccgp.fa.masked ./widow-ccgp_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/hap2/widow-ccgp.fa.masked ./widow-ccgp_2.fasta
#rohit haplotypes
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/rougheye-kolora.fa.masked ./Sebastesaleutianus_rougheye_kolora_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/rosy-kolora.fa.masked  ./Sebastesrosaceus_rosy_kolora_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/canary-kolora.fa.masked ./Sebastespinniger_canary_kolora_1.fasta
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/vermilion-kolora.fa.masked ./Sebastesminiatus_vermilion_kolora_1.fasta
#korean
cp /global/scratch/users/nicolas931010/rockfish_pangenome/repeat_annotation/repeatmasker/fasta/korean-ouc.fa.masked ./korean_1.fasta

touch assemblies.copied
```

and run `sbatch copy_assemblies.sh`

Now make a list of assemblies we can use for other scripts
```
ls *.fasta > rockfish55_assembly.list
wc -l rockfish55_assembly.list 
##55 rockfish55_assembly.list
```
