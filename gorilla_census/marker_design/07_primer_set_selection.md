Code to select a test target site for the Jan 2025 test
Katie Solari sent over a number of sets of compatible primers for us to check

Activate conda env, make directory, and upload test snp sets
```
conda activate /global/scratch/users/rdekayne/envs/geno

mkdir -p /global/scratch/users/rdekayne/gorilla_census/07_primer_set_selection && cd /global/scratch/users/rdekayne/gorilla_census/07_primer_set_selection

scp ../Katie_successful_run_output/test_set1.txt rdekayne@hpc.brc.berkeley.edu:/global/scratch/users/rdekayne/gorilla_census/07_primer_set_selection
```
Copy vcf file from mpcrselect of high pi snps and the subset to only include snps in each test set
```
cp ../05_mpcrselect/run_one_output/16_PiFinalSNPs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.* .
gunzip autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.*

##remove weird header line manually with vim

bgzip autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf
bcftools index autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz
```
Now do the subsetting and count how many SNPs each set has
```
bcftools view -R set_1_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set1_filtered.vcf
grep "chr" set1_filtered.vcf | grep -v "##" | wc -l 
#129 

bcftools view -R set_2_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set2_filtered.vcf
grep "chr" set2_filtered.vcf | grep -v "##" | wc -l 
#132 

bcftools view -R set_3_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set3_filtered.vcf
grep "chr" set3_filtered.vcf | grep -v "##" | wc -l 
#129

bcftools view -R set_4_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set4_filtered.vcf
grep "chr" set4_filtered.vcf | grep -v "##" | wc -l 
#129

bcftools view -R set_5_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set5_filtered.vcf
grep "chr" set5_filtered.vcf | grep -v "##" | wc -l 
#132

bcftools view -R set_6_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set6_filtered.vcf
grep "chr" set6_filtered.vcf | grep -v "##" | wc -l 
#129

bcftools view -R set_7_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set7_filtered.vcf
grep "chr" set7_filtered.vcf | grep -v "##" | wc -l 
#130

bcftools view -R set_8_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set8_filtered.vcf
grep "chr" set8_filtered.vcf | grep -v "##" | wc -l 
#132

bcftools view -R set_9_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set9_filtered.vcf
grep "chr" set9_filtered.vcf | grep -v "##" | wc -l 
#129

##part2
bcftools view -R part2set_1_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > part2_set1_filtered.vcf
grep "chr" part2_set1_filtered.vcf | grep -v "##" | wc -l 
#144

bcftools view -R part2set_2_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > part2_set2_filtered.vcf
grep "chr" part2_set2_filtered.vcf | grep -v "##" | wc -l 
#156
```

Now produce PCAs of each set to make sure we are differentiating individuals
```
conda activate /global/scratch/users/rdekayne/envs/pca


for file in *_filtered.vcf
do
plink --vcf ${file} --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${file}_unfilt_out
touch ./${file}.done
done

for file in part*_filtered.vcf
do
plink --vcf ${file} --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${file}_unfilt_out
touch ./${file}.done
done
```

Check overlap of different sets
```
cat set_*_loci.txt > all_sets.txt
sort -o all_sets.txt -u all_sets.txt

for file in *_filtered.vcf
do
bcftools query -f '%CHROM %POS\n' ${file} > ${file}_snp_pos.txt
done

for file in part*_filtered.vcf
do
bcftools query -f '%CHROM %POS\n' ${file} > ${file}_snp_pos.txt
done
```
```

bcftools view -R set_10_loci.txt ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.finPi.vcf.gz > set10_filtered.vcf
grep "chr" set10_filtered.vcf | grep -v "##" | wc -l 
#132
