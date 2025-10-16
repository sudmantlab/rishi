Here we will run PSMC on our assemblies

The VCF files calling variant sites has already been produced 

We will start by preparing our environment
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/06_PSMC && cd /global/scratch/users/rdekayne/rockfish55/06_PSMC

conda create -p /global/scratch/users/rdekayne/envs/psmc2
conda activate /global/scratch/users/rdekayne/envs/psmc2
conda install -c conda-forge -c bioconda bcftools=1.16
conda install conda-forge::libgcc-ng
conda install bioconda::psmc
conda install bioconda::seqtk
conda install bioconda::samtools
conda install conda-forge::gnuplot
conda install conda-forge::texlive-core
conda install conda-forge::ghostscript

conda activate /global/scratch/users/rdekayne/envs/psmc2
mkdir -p vcf_files
```

Now copy our vcf files here and remove the black and yellow rockfish which was not analysed due to being poor quality
```
cp /global/scratch/users/nicolas931010/rockfish_pangenome/sv_detection/paftools/hap2_to_hap1/*.vcf ./vcf_files/.
rm ./vcf_files/SebasteschrysomelasSEB-260.vcf
ls ./vcf_files/* > indiv_list.txt
```

This is what our `indiv_list.txt` looks like
```
./vcf_files/7164-PS-0001.vcf
./vcf_files/7164-PS-0002.vcf
./vcf_files/7164-PS-0003.vcf
./vcf_files/7164-PS-0004.vcf
./vcf_files/7164-PS-0005.vcf
./vcf_files/7164-PS-0006.vcf
./vcf_files/7164-PS-0007.vcf
./vcf_files/7164-PS-0008.vcf
./vcf_files/acadian-cbp.vcf
./vcf_files/acadian-vgp.vcf
./vcf_files/bocaccio-ccgp.vcf
./vcf_files/honeycomb-vgp.vcf
./vcf_files/m84185_240412_222911_s4.vcf
./vcf_files/SebastesalascanusSEB-2.vcf
./vcf_files/SebastesalutusPOP6.vcf
./vcf_files/SebastesborealisSEB-8.vcf
./vcf_files/SebastescarnatusSEB-258.vcf
./vcf_files/SebastescrocotulusSEB-252.vcf
./vcf_files/SebastesmelanostictusSEB-4.vcf
./vcf_files/SebastesminiatusSEB-256.vcf
./vcf_files/SebastesmystinusSEB-254.vcf
./vcf_files/SebastespolyspinisSEB-5.vcf
./vcf_files/Sebastesruberrimus12-Yellowye1a.vcf
./vcf_files/SebastesserranoidesSEB-259.vcf
./vcf_files/widow-ccgp.vcf
```

Now make some sed changes to change the filepaths
```
sed -i 's/\.\/vcf_files\///g' indiv_list.txt
sed -i 's/\.vcf//g' indiv_list.txt

#now get a list of the locations of the same .fa.masked assemblies
ls -lh ../assemblies/*_1.fasta
```

This is what our `fasta_list.txt` looks like
```
7164-PS-0001_copper_1.fasta
7164-PS-0002_black_1.fasta
7164-PS-0003_rougheye_1.fasta
7164-PS-0004_widow_1.fasta
7164-PS-0005_yellowtail_1.fasta
7164-PS-0006_popb_1.fasta
7164-PS-0007_lightdusky_1.fasta
7164-PS-0008_quillback_1.fasta
acadian-cbp_1.fasta
acadian-vgp_1.fasta
bocaccio-ccgp_1.fasta
honeycomb-vgp_1.fasta
m84185_240412_222911_s4_deacon_1.fasta
SebastesalascanusSEB-2_shortspinethornyhead_1.fasta
SebastesalutusPOP6_popa_1.fasta
SebastesborealisSEB-8_shortraker_1.fasta
SebastescarnatusSEB-258_gopher_1.fasta
SebastescrocotulusSEB-252_sunset_1.fasta
SebastesmelanostictusSEB-4_blackspotted_1.fasta
SebastesminiatusSEB-256_vermilion_1.fasta
SebastesmystinusSEB-254_blue_1.fasta
SebastespolyspinisSEB-5_northern_1.fasta
Sebastesruberrimus12-Yellowye1a_1.fasta
SebastesserranoidesSEB-259_olive_1.fasta
widow-ccgp_1.fasta
```

Now we will run all assemblies through psmc using the same paramters using `all_psmc.sh` - later we will plot with species-specific values but will use the same .psmc output which are the key files
```
#!/bin/bash
#SBATCH --job-name=psmc
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio4_htc
#SBATCH --output=psmc%j.out # output file
#SBATCH --error=psmc%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=24 # 2 CPUs per job

numb=${SLURM_ARRAY_TASK_ID}
file_name=$(cat /global/scratch/users/rdekayne/rockfish55/06_PSMC/indiv_list.txt | sed -n ${numb}p)
fasta_name=$(cat /global/scratch/users/rdekayne/rockfish55/06_PSMC/fasta_list.txt | sed -n ${numb}p)

sed 's/1\/1/0\/1/g' /global/scratch/users/nicolas931010/rockfish_pangenome/sv_detection/paftools/hap2_to_hap1/${file_name}.vcf | bcftools +fill-tags -o ./${file_name}.filt.vcf -- -t AN,AC,AF 
bcftools view -m2 -M2 -v snps ${file_name}.filt.vcf | /global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/vcf2snp.pl - > ${file_name}.snp
seqtk mutfa ../../assemblies/${fasta_name} ${file_name}.snp > ${file_name}.pseudo.fa
seqtk seq -cM /global/scratch/users/nicolas931010/rockfish_pangenome/sv_detection/paftools/hap2_to_hap1/${file_name}.callable.bed -l80 ${file_name}.pseudo.fa > ${file_name}.masked.fa
seqtk seq -l0 ${file_name}.masked.fa > ${file_name}.masked.seq.fa
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/fq2psmcfa ${file_name}.masked.seq.fa > ${file_name}.psmcfa
seqtk seq -l0 ${file_name}.psmcfa > ${file_name}.seq.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${file_name}.psmc ${file_name}.seq.psmcfa
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.36e-09 -g 7.4 -p ${file_name}.plot ${file_name}.psmc

touch ${file_name}.done
```
And now run this
```
wc -l indiv_list.txt 
##25 indiv_list.txt
sbatch --array=1-25 all_psmc.sh
```

Now we are going to plot our psmc output files
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/06_PSMC/combined 
cp *.psmc /global/scratch/users/rdekayne/rockfish55/06_PSMC/combined/
cd /global/scratch/users/rdekayne/rockfish55/06_PSMC/combined
rm honeycomb-vgp.psmc

cat *.psmc > combined.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -M '7164-PS-0001.psmc,7164-PS-0002.psmc,7164-PS-0003.psmc,7164-PS-0004.psmc,7164-PS-0005.psmc,7164-PS-0006.psmc,7164-PS-0007.psmc,7164-PS-0008.psmc,acadian-cbp.psmc,acadian-vgp.psmc,bocaccio-ccgp.psmc,m84185_240412_222911_s4.psmc,SebastesalascanusSEB-2.psmc,SebastesalutusPOP6.psmc,SebastesborealisSEB-8.psmc,SebastescarnatusSEB-258.psmc,SebastescrocotulusSEB-252.psmc,SebastesmelanostictusSEB-4.psmc,SebastesminiatusSEB-256.psmc,SebastesmystinusSEB-254.psmc,SebastespolyspinisSEB-5.psmc,Sebastesruberrimus12-Yellowye1a.psmc,SebastesserranoidesSEB-259.psmc,widow-ccgp.psmc' -u 4.36e-09 -g 7.4 -f Helvetica,8 -p combined_plot.plot combined.psmc
```

Now do this again but I'm going to modify the plotting script to minimise legend font size
```
cp /global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl ./psmc_plot_UPDATED.pl
my $legend_font_size = 8; # Set this to your desired legend font size
my $keyconf = $opts{M} ? "set key $opts{P} font 'Helvetica,$legend_font_size'" : "set key off";
```

And now run this updated script
```
chmod +x psmc_plot_UPDATED.pl
./psmc_plot_UPDATED.pl -M 'copper,black,rougheye,widow,yellowtail,pop-B,light-dusky,quillback,acadian-cbp,acadian-vgp,bocaccio,honeycomb-vgp,deacon,shortspine-thornyhead,pop-A,shortraker,gopher,sunset,blackspotted,vermilion,blue,northern,yelloweye,olive,widow' -u 4.36e-09 -g 7.4 -Y 45 -p combined_plot_UPDATED.plot combined.psmc
```

Now we will produce one plot per main clade in our tree

Prepare the environment
```
conda activate /global/scratch/users/rdekayne/envs/psmc2
mkdir -p /global/scratch/users/rdekayne/rockfish55/06_PSMC/clade_psmc && cd /global/scratch/users/rdekayne/rockfish55/06_PSMC/clade_psmc
cp ../combined/psmc_plot_UPDATED.pl .
```

The first clade 'GROUP1' will include
- seb8 shortraker
- seb4 blackspotted
- 0003 rougheye
- seb5 northern
- 0007 lightdusky
- pop6 popa
- 0006 popb
- acadian

So we combine the psmc files and plot
```
cat ../SebastesborealisSEB-8.psmc ../SebastesmelanostictusSEB-4.psmc ../7164-PS-0003.psmc ../SebastespolyspinisSEB-5.psmc ../7164-PS-0007.psmc ../SebastesalutusPOP6.psmc ../7164-PS-0006.psmc ../acadian-cbp.psmc ../acadian-vgp.psmc > group1.psmc

./psmc_plot_UPDATED.pl -M 'shortraker,blackspotted,rougheye,northern,lightdusky,popA,popB,acadian-cbp,acadian-vgp' -u 4.36e-09 -g 22.7 -Y 45 -p group1_22.7.plot group1.psmc
```

The second clade 'GROUP2' includes
- seb259 olive
- 0002 black
- 0005 yellowtail
- widow
- 0004 widow
- seb254 blue
- m84 deacon

And again cat and plot
```
cat ../SebastesserranoidesSEB-259.psmc ../7164-PS-0002.psmc ../7164-PS-0005.psmc ../widow-ccgp.psmc ../7164-PS-0004.psmc ../SebastesmystinusSEB-254.psmc ../m84185_240412_222911_s4.psmc > group2.psmc

./psmc_plot_UPDATED.pl -M 'olive,black,yellowtail,widow-ccgp,widow,blue,deacon' -u 4.36e-09 -g 22.7 -Y 45 -p group2_22.7.plot group2.psmc
```

And the final group 'GROUP3' includes
- seb252 sunset
- seb256 vermilion
- bocaccio
- yelloweye
- honeycomb
- seb258 gopher
- 0001 copper
- 0008 quillback

```
cat ../SebastescrocotulusSEB-252.psmc ../SebastesminiatusSEB-256.psmc ../bocaccio-ccgp.psmc ../Sebastesruberrimus12-Yellowye1a.psmc ../honeycomb-vgp.psmc ../SebastescarnatusSEB-258.psmc ../7164-PS-0001.psmc ../7164-PS-0008.psmc > group3.psmc

./psmc_plot_UPDATED.pl -M 'sunset,vermilion,bocaccio-ccgp,yelloweye,honeycomb-vgp,gopher,copper,quillback' -u 4.36e-09 -g 22.7 -Y 45 -p group3_22.7.plot group3.psmc
```

Now plot them all on the same scale
```
./psmc_plot_UPDATED.pl -M 'shortraker,blackspotted,rougheye,northern,lightdusky,popA,popB,acadian-cbp,acadian-vgp' -u 4.36e-09 -g 22.7 -Y 45 -X 10000000 -p group1_22.7_scaled.plot group1.psmc
./psmc_plot_UPDATED.pl -M 'olive,black,yellowtail,widow-ccgp,widow,blue,deacon' -u 4.36e-09 -g 22.7 -Y 45 -X 10000000 -p group2_22.7_scaled.plot group2.psmc
./psmc_plot_UPDATED.pl -M 'sunset,vermilion,bocaccio-ccgp,yelloweye,honeycomb-vgp,gopher,copper,quillback' -u 4.36e-09 -g 22.7 -Y 45 -X 10000000 -p group3_22.7_scaled.plot group3.psmc
```

Now we will run our generation-time-specific plots (these will then be combined in an R script)
```
mkdir -p /global/scratch/users/rdekayne/rockfish55/06_PSMCgen_times && cd /global/scratch/users/rdekayne/rockfish55/06_PSMC/gen_times

###G1
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 44.89 -p SebastesborealisSEB-8_shortraker.plot ../SebastesborealisSEB-8.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 44.89 -p SebastesmelanostictusSEB-4_blackspotted.plot ../SebastesmelanostictusSEB-4.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 44.89 -p 7164-PS-0003_rougheye.plot ../7164-PS-0003.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p SebastespolyspinisSEB-5_northern.plot ../SebastespolyspinisSEB-5.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p 7164-PS-0007_lightdusky.plot ../7164-PS-0007.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p SebastesalutusPOP6_popA.plot ../SebastesalutusPOP6.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p 7164-PS-0006_popB.plot ../7164-PS-0006.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p acadian-cbp.plot ../acadian-cbp.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 21.55 -p acadian-vgp.plot ../acadian-vgp.psmc

###G2
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 9.38 -p SebastesserranoidesSEB-259_olive.plot ../SebastesserranoidesSEB-259.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.69 -p 7164-PS-0002_black.plot ../7164-PS-0002.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 16.99 -p 7164-PS-0005_yellowtail.plot ../7164-PS-0005.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.37 -p widow-ccgp.plot ../widow-ccgp.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.37 -p 7164-PS-0004_widow.plot ../7164-PS-0004.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.38 -p SebastesmystinusSEB-254_blue.plot ../SebastesmystinusSEB-254.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.38 -p m84185_240412_222911_s4_deacon.plot ../m84185_240412_222911_s4.psmc

###G3
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 20.20 -p SebastescrocotulusSEB-252_sunset.plot ../SebastescrocotulusSEB-252.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 20.20 -p SebastesminiatusSEB-256_vermilion.plot ../SebastesminiatusSEB-256.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 20.20 -p bocaccio-ccgp.plot ../bocaccio-ccgp.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 40.72 -p Sebastesruberrimus12-Yellowye1a.plot ../Sebastesruberrimus12-Yellowye1a.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 20.20 -p honeycomb-vgp.plot ../honeycomb-vgp.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 20.20 -p SebastescarnatusSEB-258_gopher.plot ../SebastescarnatusSEB-258.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 14.97 -p 7164-PS-0001_copper.plot ../7164-PS-0001.psmc
/global/scratch/users/rdekayne/rockfish55/06_PSMC/psmc/utils/psmc_plot.pl -R -u 4.743333e-09 -g 25.42 -p 7164-PS-0008_quillback.plot ../7164-PS-0008.psmc
```

Then in R we run `plot_psmc.R` as follows
```
#script to plot PSMC output for rockfish - De-Kayne 2025

#remotes::install_github("emmanuelparadis/psmcr/psmcr")
library(psmcr)
library(tidyverse)
library(cols4all)

setwd('/Users/rdk_mac_work/Dropbox/RishiMac3Work/rockfish/rockfish55/06_PSMC/')


# Read in PSMC results ----------------------------------------------------

##G1
acadiancbp <- import.psmc('acadian-cbp.psmc')
acadianvgp <- import.psmc('acadian-vgp.psmc')
rougheye <- import.psmc('7164-PS-0003.psmc')
popb <- import.psmc('7164-PS-0006.psmc')
lightdusky <- import.psmc('7164-PS-0007.psmc')
popa <- import.psmc('SebastesalutusPOP6.psmc')
shortraker <- import.psmc('SebastesborealisSEB-8.psmc')
blackspotted <- import.psmc('SebastesmelanostictusSEB-4.psmc')
northern <- import.psmc('SebastespolyspinisSEB-5.psmc')

##G2
black <- import.psmc('7164-PS-0002.psmc')
widow <- import.psmc('7164-PS-0004.psmc')
yellowtail <- import.psmc('7164-PS-0005.psmc')
blue <- import.psmc('SebastesmystinusSEB-254.psmc')
olive <- import.psmc('SebastesserranoidesSEB-259.psmc')
deacon <- import.psmc('m84185_240412_222911_s4.psmc')
widowccgp <- import.psmc('widow-ccgp.psmc')


##G3
copper<- import.psmc('7164-PS-0001.psmc')
quillback <- import.psmc('7164-PS-0008.psmc')
gopher <- import.psmc('SebastescarnatusSEB-258.psmc')
sunset <- import.psmc('SebastescrocotulusSEB-252.psmc')
vermilion <- import.psmc('SebastesminiatusSEB-256.psmc')
yelloweye <- import.psmc('Sebastesruberrimus12-Yellowye1a.psmc')
bocaccioccgp <- import.psmc('bocaccio-ccgp.psmc')
honeycombvgp <- import.psmc('honeycomb-vgp.psmc')


#def variables
#mutation.rate <- 4.56e-9
#g  <-  1               # generation set to 1, already taken into account in mutation rate value
scaled <- F
bin.size <- 100
ylim <- NULL
xlim <- NULL


#functions to scale and extract psmc data
.getXYplot.psmc <- function(x, mutation.rate, g, scaled, bin.size){
  RS <- x$RS[x$RS[, "iter"] == x$niters, ]
  theta0 <- x$parameters[nrow(x$parameters), "theta0"]
  if (scaled) {
    ##xx <- RS[, "t_k"] / (theta0 / bin.size)
    xx <- RS[, "t_k"] / (theta0 * bin.size)
    yy <- theta0 * RS[, "lambda_k"] / bin.size
  } else {
    ##N0 <- theta0/(4 * mutation.rate / bin.size)
    denom <- 4 * mutation.rate * g * bin.size
    N0 <- theta0 / denom
    xx <- 2 * N0 * RS[, "t_k"]
    yy <- N0 * RS[, "lambda_k"]
  }
  list(xx = xx, yy = yy)
}

getPsmcPlotData <- function(x, species.name, mutation.rate, g, scaled=F, bin.size, xlim=NULL, ylim=NULL) {
  xy <- .getXYplot.psmc(x, mutation.rate = mutation.rate, g = g,
                        scaled = scaled, bin.size = bin.size)
  xx <- xy$xx
  yy <- xy$yy
  
  if (is.null(ylim)) yl <- max(yy)
  withbootstrap <- !is.null(x$bootstrap)
  if (withbootstrap) {
    obj <- .getBootstrap.psmc(x, mutation.rate = mutation.rate, g = g,
                              scaled = scaled, bin.size = bin.size)
    Tk.boot <- obj$Tk.boot
    Nk.boot <- obj$Nk.boot
    if (is.null(ylim)) yl <- max(yl, Nk.boot)
  }
  if (is.null(xlim)) xlim <- range(xx)
  if (is.null(ylim)) ylim <- c(0, yl)
  
  outdata <- data.frame('x'=xx,'y'=yy,'species'=species.name)
}



#get plot data for each species
#G1
acadiancbp.plot <- getPsmcPlotData(acadiancbp,'acadian-cbp',4.36e-9,21.55,scaled,bin.size,ylim,xlim)
acadianvgp.plot <- getPsmcPlotData(acadianvgp,'acadian-vgp',4.36e-9,21.55,scaled,bin.size,ylim,xlim)
rougheye.plot <- getPsmcPlotData(rougheye,'rougheye',4.36e-9,44.89,scaled,bin.size,ylim,xlim)
popa.plot <- getPsmcPlotData(popa,'POP A',4.36e-9,21.55,scaled,bin.size,ylim,xlim)
popb.plot <- getPsmcPlotData(popb,'POP B',4.36e-9,21.55,scaled,bin.size,ylim,xlim)
lightdusky.plot <- getPsmcPlotData(lightdusky,'lightdusky',4.36e-9,21.55,scaled,bin.size,ylim,xlim)
shortraker.plot <- getPsmcPlotData(shortraker,'shortraker',4.36e-9,44.89,scaled,bin.size,ylim,xlim)
blackspotted.plot <- getPsmcPlotData(blackspotted,'blackspotted',4.36e-9,44.89,scaled,bin.size,ylim,xlim)
northern.plot <- getPsmcPlotData(northern,'northern',4.36e-9,44.89,scaled,bin.size,ylim,xlim)

#G2
black.plot <- getPsmcPlotData(black,'black',4.36e-9,14.69,scaled,bin.size,ylim,xlim)
widow.plot <- getPsmcPlotData(widow,'widow',4.36e-9,14.37,scaled,bin.size,ylim,xlim)
yellowtail.plot <- getPsmcPlotData(yellowtail,'yellowtail',4.36e-9,16.99,scaled,bin.size,ylim,xlim)
blue.plot <- getPsmcPlotData(blue,'blue',4.36e-9,14.38,scaled,bin.size,ylim,xlim)
olive.plot <- getPsmcPlotData(olive,'olive',4.36e-9,9.38,scaled,bin.size,ylim,xlim)
deacon.plot <- getPsmcPlotData(deacon,'deacon',4.36e-9,14.38,scaled,bin.size,ylim,xlim)
widowccgp.plot <- getPsmcPlotData(widowccgp,'widowccgp',4.36e-9,14.37,scaled,bin.size,ylim,xlim)

#G3
copper.plot <- getPsmcPlotData(copper,'copper',4.36e-9,14.97,scaled,bin.size,ylim,xlim)
quillback.plot <- getPsmcPlotData(quillback,'quillback',4.36e-9,25.42,scaled,bin.size,ylim,xlim)
gopher.plot <- getPsmcPlotData(gopher,'gopher',4.36e-9,20.20,scaled,bin.size,ylim,xlim)
sunset.plot <- getPsmcPlotData(sunset,'sunset',4.36e-9,20.20,scaled,bin.size,ylim,xlim)
vermilion.plot <- getPsmcPlotData(vermilion,'vermilion',4.36e-9,20.20,scaled,bin.size,ylim,xlim)
yelloweye.plot <- getPsmcPlotData(yelloweye,'yelloweye',4.36e-9,40.72,scaled,bin.size,ylim,xlim)
bocaccioccgp.plot <- getPsmcPlotData(bocaccioccgp,'bocaccioccgp',4.36e-9,20.20,scaled,bin.size,ylim,xlim)
honeycombvgp.plot <- getPsmcPlotData(honeycombvgp,'honeycombvgp',4.36e-9,20.20,scaled,bin.size,ylim,xlim)

#####
acadiancbp.plot <- getPsmcPlotData(acadiancbp,'acadian-cbp',4.36e-9,1,scaled,bin.size,ylim,xlim)
acadianvgp.plot <- getPsmcPlotData(acadianvgp,'acadian-vgp',4.36e-9,1,scaled,bin.size,ylim,xlim)
rougheye.plot <- getPsmcPlotData(rougheye,'rougheye',4.36e-9,1,scaled,bin.size,ylim,xlim)
popa.plot <- getPsmcPlotData(popa,'POP A',4.36e-9,1,scaled,bin.size,ylim,xlim)
popb.plot <- getPsmcPlotData(popb,'POP B',4.36e-9,1,scaled,bin.size,ylim,xlim)
lightdusky.plot <- getPsmcPlotData(lightdusky,'lightdusky',4.36e-9,1,scaled,bin.size,ylim,xlim)
shortraker.plot <- getPsmcPlotData(shortraker,'shortraker',4.36e-9,1,scaled,bin.size,ylim,xlim)
blackspotted.plot <- getPsmcPlotData(blackspotted,'blackspotted',4.36e-9,1,scaled,bin.size,ylim,xlim)
northern.plot <- getPsmcPlotData(northern,'northern',4.36e-9,1,scaled,bin.size,ylim,xlim)

#G2
black.plot <- getPsmcPlotData(black,'black',4.36e-9,1,scaled,bin.size,ylim,xlim)
widow.plot <- getPsmcPlotData(widow,'widow',4.36e-9,1,scaled,bin.size,ylim,xlim)
yellowtail.plot <- getPsmcPlotData(yellowtail,'yellowtail',4.36e-9,1,scaled,bin.size,ylim,xlim)
blue.plot <- getPsmcPlotData(blue,'blue',4.36e-9,1,scaled,bin.size,ylim,xlim)
olive.plot <- getPsmcPlotData(olive,'olive',4.36e-9,1,scaled,bin.size,ylim,xlim)
deacon.plot <- getPsmcPlotData(deacon,'deacon',4.36e-9,1,scaled,bin.size,ylim,xlim)
widowccgp.plot <- getPsmcPlotData(widowccgp,'widowccgp',4.36e-9,1,scaled,bin.size,ylim,xlim)

#G3
copper.plot <- getPsmcPlotData(copper,'copper',4.36e-9,1,scaled,bin.size,ylim,xlim)
quillback.plot <- getPsmcPlotData(quillback,'quillback',4.36e-9,1,scaled,bin.size,ylim,xlim)
gopher.plot <- getPsmcPlotData(gopher,'gopher',4.36e-9,1,scaled,bin.size,ylim,xlim)
sunset.plot <- getPsmcPlotData(sunset,'sunset',4.36e-9,1,scaled,bin.size,ylim,xlim)
vermilion.plot <- getPsmcPlotData(vermilion,'vermilion',4.36e-9,1,scaled,bin.size,ylim,xlim)
yelloweye.plot <- getPsmcPlotData(yelloweye,'yelloweye',4.36e-9,1,scaled,bin.size,ylim,xlim)
bocaccioccgp.plot <- getPsmcPlotData(bocaccioccgp,'bocaccioccgp',4.36e-9,1,scaled,bin.size,ylim,xlim)
honeycombvgp.plot <- getPsmcPlotData(honeycombvgp,'honeycombvgp',4.36e-9,1,scaled,bin.size,ylim,xlim)


# Combine all species into single dataframe for multi-plot ----------------
##G1
combined.plot <- acadiancbp.plot %>% 
  bind_rows(acadianvgp.plot, rougheye.plot, popa.plot, popb.plot, lightdusky.plot, shortraker.plot, blackspotted.plot, northern.plot)

G1_PLOT <- ggplot(combined.plot,aes(x=x+1,y=y/1e4,color=species)) +
            geom_step(lwd=1) +
            scale_x_log10() +
            coord_cartesian(xlim = c(1.5e4,1.5e9), ylim = c(0,10000)) +
            scale_color_discrete_c4a_cat(palette = 'safe') +
            scale_fill_brewer(palette = 'Set2') +
            labs(x='Time',y='N (x10^4)',color='Species') +
            theme_linedraw() + theme(legend.text = element_text(face='italic'),strip.background = element_rect(fill='grey10'),
                                     strip.text = element_text(face='italic'),
                                     panel.grid = element_blank())
G1_PLOT
ggsave("G1.pdf", plot = G1_PLOT, width = 10, height = 6) 


##G2
combined.plot <- black.plot %>% 
  bind_rows(widow.plot, yellowtail.plot, blue.plot, olive.plot, deacon.plot, widowccgp.plot)

G2_PLOT <- ggplot(combined.plot,aes(x=x+1,y=y/1e4,color=species)) +
  geom_step(lwd=1) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1.5e1,1.5e5), ylim = c(0,1.5)) +
  scale_color_discrete_c4a_cat(palette = 'safe') +
  scale_fill_brewer(palette = 'Set2') +
  labs(x='Time',y='N (x10^4)',color='Species') +
  theme_linedraw() + theme(legend.text = element_text(face='italic'),strip.background = element_rect(fill='grey10'),
                           strip.text = element_text(face='italic'),
                           panel.grid = element_blank())
G2_PLOT
ggsave("G2.pdf", plot = G2_PLOT, width = 10, height = 6) 

##G3
combined.plot <- copper.plot %>% 
  bind_rows(quillback.plot, gopher.plot, sunset.plot, vermilion.plot, yelloweye.plot, bocaccioccgp.plot, honeycombvgp.plot)

G3_PLOT <- ggplot(combined.plot,aes(x=x+1,y=y/1e4,color=species)) +
  geom_step(lwd=1) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1.5e1,1.5e5), ylim = c(0,1.5)) +
  scale_color_discrete_c4a_cat(palette = 'safe') +
  scale_fill_brewer(palette = 'Set2') +
  labs(x='Time',y='N (x10^4)',color='Species') +
  theme_linedraw() + theme(legend.text = element_text(face='italic'),strip.background = element_rect(fill='grey10'),
                           strip.text = element_text(face='italic'),
                           panel.grid = element_blank())
G3_PLOT
ggsave("G3.pdf", plot = G3_PLOT, width = 10, height = 6) 

##ALL
combined.plot <- copper.plot %>% 
  bind_rows(quillback.plot, gopher.plot, sunset.plot, vermilion.plot, yelloweye.plot, bocaccioccgp.plot, honeycombvgp.plot,
            acadiancbp.plot, acadianvgp.plot, rougheye.plot, popa.plot, popb.plot, lightdusky.plot, shortraker.plot, blackspotted.plot, northern.plot,
            black.plot, widow.plot, yellowtail.plot, blue.plot, olive.plot, deacon.plot, widowccgp.plot)

ggplot(combined.plot,aes(x=x+1,y=y/1e4,color=species)) +
  geom_step(lwd=1) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1.5e1,1.5e5), ylim = c(0,1.5)) +
  #scale_color_discrete_c4a_cat(palette = 'safe') +
  scale_fill_brewer(palette = 'Set2') +
  labs(x='Time',y='N (x10^4)',color='Species') +
  theme_linedraw() + theme(legend.text = element_text(face='italic'),strip.background = element_rect(fill='grey10'),
                           strip.text = element_text(face='italic'),
                           panel.grid = element_blank())


#######################################################################
#single species plot

ggplot(acadiancbp.plot,aes(x=x+1,y=y/1e4,color=species)) +
  geom_step(lwd=1) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1.5e1,1.5e5), ylim = c(0,1.5)) +
  # scale_color_discrete_c4a_cat(palette = 'safe') +
  scale_fill_brewer(palette = 'Set2') +
  labs(x='Time',y='N (x10^4)',color='Species') +
  theme_linedraw() + theme(legend.text = element_text(face='italic'),strip.background = element_rect(fill='grey10'),
                           strip.text = element_text(face='italic'),
                           panel.grid = element_blank())
```

![psmcplots](plots/psmc_plots.pdf "psmc plots")
