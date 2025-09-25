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
```

Decided to select part2set_2_loci.txt which had 156 loci - `part2set_2_loci.txt `
```
chr1pathsa1	16135939	16135988
chr1pathsa1	21927278	21927350
chr1pathsa1	30688507	30688600
chr1pathsa1	38452429	38452503
chr1pathsa1	54998354	54998419
chr1pathsa1	63234873	63234944
chr1pathsa1	64133967	64134057
chr1pathsa1	88152097	88152155
chr1pathsa1	89063324	89063373
chr1pathsa1	91602073	91602175
chr1pathsa1	119974287	119974343
chr1pathsa1	120782994	120783043
chr1pathsa1	136905386	136905452
chr1pathsa1	150725111	150725165
chr1pathsa1	157979018	157979108
chr1pathsa1	160543318	160543403
chr1pathsa1	161048913	161048975
chr1pathsa1	172049091	172049140
chr1pathsa1	176066464	176066532
chr1pathsa1	202329613	202329662
chr1pathsa1	234606339	234606405
chr2pathsa3	10341626	10341677
chr2pathsa3	15697826	15697894
chr2pathsa3	18333120	18333199
chr2pathsa3	19191865	19191928
chr2pathsa3	36437375	36437431
chr2pathsa3	44368885	44368943
chr2pathsa3	66182308	66182376
chr2pathsa3	72201535	72201617
chr2pathsa3	73618155	73618212
chr2pathsa3	88314689	88314746
chr2pathsa3	92599726	92599791
chr2pathsa3	107300409	107300482
chr2pathsa3	122253260	122253332
chr2pathsa3	139497068	139497134
chr2pathsa3	146849051	146849100
chr2pathsa3	158501194	158501246
chr2pathsa3	188644409	188644484
chr2pathsa3	190070017	190070067
chr3pathsa4	25746794	25746846
chr3pathsa4	28278972	28279029
chr3pathsa4	41166231	41166280
chr3pathsa4	90175751	90175803
chr3pathsa4	100143789	100143839
chr3pathsa4	106849390	106849440
chr3pathsa4	112752270	112752333
chr3pathsa4	115422765	115422815
chr3pathsa4	129489391	129489491
chr3pathsa4	169046587	169046636
chr3pathsa4	176784175	176784282
chr3pathsa4	199792180	199792229
chr3pathsa4	210047762	210047811
chr3pathsa4	210350525	210350600
chr4pathsa17x5	17872791	17872841
chr4pathsa17x5	25667135	25667208
chr4pathsa17x5	26072463	26072512
chr4pathsa17x5	35934900	35934950
chr4pathsa17x5	57210711	57210765
chr4pathsa17x5	96726263	96726313
chr4pathsa17x5	99427365	99427414
chr4pathsa17x5	105669207	105669300
chr4pathsa17x5	118525887	118525950
chr4pathsa17x5	126606594	126606648
chr4pathsa17x5	131579671	131579745
chr4pathsa17x5	134691587	134691636
chr4pathsa17x5	144731517	144731566
chr4pathsa17x5	155211167	155211232
chr4pathsa17x5	172000676	172000728
chr5mathsa6	32015973	32016049
chr5mathsa6	43593146	43593199
chr5mathsa6	59190327	59190381
chr5mathsa6	65503372	65503421
chr5mathsa6	71307165	71307219
chr5mathsa6	93747683	93747752
chr5mathsa6	95573766	95573823
chr5mathsa6	113852937	113853012
chr5mathsa6	117155646	117155707
chr5mathsa6	120423745	120423804
chr5mathsa6	137881526	137881636
chr5mathsa6	157062069	157062134
chr5mathsa6	185994009	185994061
chr6mathsa7	38658870	38658927
chr6mathsa7	91703053	91703102
chr6mathsa7	91804192	91804242
chr6mathsa7	98938729	98938778
chr6mathsa7	108809554	108809671
chr6mathsa7	121100608	121100664
chr6mathsa7	149030557	149030633
chr6mathsa7	155486551	155486617
chr6mathsa7	168040177	168040247
chr7pathsa8	8880800	8880853
chr7pathsa8	20386282	20386344
chr7pathsa8	39190975	39191024
chr7pathsa8	45067078	45067127
chr7pathsa8	53147698	53147758
chr7pathsa8	66606033	66606081
chr7pathsa8	86979745	86979798
chr7pathsa8	89270764	89270820
chr7pathsa8	103059304	103059363
chr7pathsa8	105275710	105275759
chr7pathsa8	107897717	107897766
chr7pathsa8	127041663	127041722
chr8pathsa10	38948573	38948630
chr8pathsa10	92940203	92940254
chr8pathsa10	108955326	108955378
chr8pathsa10	120898123	120898202
chr8pathsa10	128120986	128121033
chr8pathsa10	138671457	138671507
chr9pathsa11	18434776	18434835
chr9pathsa11	32322872	32322948
chr9pathsa11	41286296	41286386
chr9pathsa11	44312794	44312850
chr9pathsa11	48976736	48976799
chr9pathsa11	50101421	50101470
chr9pathsa11	86688494	86688540
chr9pathsa11	87190721	87190771
chr9pathsa11	88212087	88212140
chr9pathsa11	90673597	90673660
chr9pathsa11	92087786	92087835
chr9pathsa11	108474861	108474941
chr9pathsa11	109498618	109498678
chr9pathsa11	110525519	110525568
chr9pathsa11	136877075	136877138
chr10mathsa12	72338598	72338648
chr11mathsa2b	69948297	69948350
chr11mathsa2b	76168504	76168573
chr11mathsa2b	87171702	87171809
chr11mathsa2b	112907540	112907632
chr11mathsa2b	122862941	122862995
chr11mathsa2b	123991767	123991839
chr12pathsa2a	64419143	64419214
chr12pathsa2a	67368132	67368197
chr12pathsa2a	97933263	97933314
chr12pathsa2a	129058900	129058974
chr13pathsa9	68582754	68582843
chr13pathsa9	74812423	74812489
chr13pathsa9	75416496	75416580
chr13pathsa9	113520582	113520652
chr13pathsa9	118847672	118847728
chr14pathsa13	99776586	99776650
chr15pathsa14	123256877	123256953
chr16pathsa15	55895902	55895952
chr17mathsa18	41800590	41800660
chr17mathsa18	84483023	84483078
chr17mathsa18	85290203	85290294
chr17mathsa18	103963825	103963877
chr18pathsa16	81043713	81043789
chr18pathsa16	85131689	85131756
chr18pathsa16	95353497	95353551
chr18pathsa16	102360935	102360984
chr19pathsa5x17	99139724	99139781
chr19pathsa5x17	105203101	105203156
chr21pathsa20	47436332	47436412
chr21pathsa20	73112857	73112912
chr22mathsa21	43262977	43263026
chr23mathsa22	46384383	46384433
```

