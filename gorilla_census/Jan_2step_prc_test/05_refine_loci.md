Code for selecting new loci for second test using a 1step pcr to reduce risk of contamination

See how many loci from 40 selected SNPs (pickeda at random) we sequenced here

Target amplicons in: `40_singlestep_amplicons.txt`

Calculate counts overlapping using our target amplicons as a regions file
```
bgzip targets156_specific_filtered_MAC0.vcf
bcftools index targets156_specific_filtered_MAC0.vcf.gz
bcftools view -R 40_singlestep_amplicons.txt targets156_specific_filtered_MAC0.vcf.gz > targets156_specific_filtered_MAC0_40ss_amplicons.vcf
bcftools view -H targets156_specific_filtered_MAC0_40ss_amplicons.vcf | wc -l
```
32

Check statistics for the kit we selected (low volume)
```
cp sml_vol_sample_list.txt sml_vol_sample_list_NO6.txt
#remove well6 == negative control
bcftools view -S sml_vol_sample_list_NO6.txt -o targets156_specific_filtered_MAC0_40ss_amplicons_SML_no6.vcf targets156_specific_filtered_MAC0_40ss_amplicons.vcf

vcftools --vcf targets156_specific_filtered_MAC0_40ss_amplicons_SML_no6.vcf --missing-site --out SML_40ss_no6_
```

Now identify high pi SNPs to send to Katie Solari for new loci creation to get 40 high diversity SNPs and confirm the SNPs will not amplify in humans
```
mkdir -p /global/scratch/users/rdekayne/gorilla_census/12_human_primer_check && cd /global/scratch/users/rdekayne/gorilla_census/12_human_primer_check

cp /global/scratch/users/rdekayne/gorilla_census/05_mpcrselect/run_one_output/10_PlinkLD/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.pruned.vcf.gz .

vcftools --gzvcf ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.pruned.vcf.gz --site-pi --out pi_output

awk 'NR==1 || $3 >= 0.25' pi_output.sites.pi > pi_output.sites.pi_over_25.pi

cut -f1,2 pi_output.sites.pi_over_25.pi > positions.txt

vcftools --gzvcf ./autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.pruned.vcf.gz --positions positions.txt --recode --recode-INFO-all --out gorilla_snps_with_pi_over_25.vcf
```
And share this file
