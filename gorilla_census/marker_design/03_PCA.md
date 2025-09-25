Make PCAs from combined vcf

Make directory and conda environments
```
conda create -p /global/scratch/users/rdekayne/envs/pca
conda activate /global/scratch/users/rdekayne/envs/pca
conda install bioconda::plink

mkdir -p /global/scratch/users/rdekayne/gorilla_census/03_PCA && cd /global/scratch/users/rdekayne/gorilla_census/03_PCA
```
Now run plink to create PCAs which we plot later R - `03.1_PCA_unfilt.sh`
```
#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --time=0-4:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3
#SBATCH --output=pca.%j.out # output file
#SBATCH --error=pca.%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=1 # 2 CPUs per job

cd /global/scratch/users/rdekayne/gorilla_census/03_PCA

plink --vcf /global/scratch/users/rdekayne/gorilla_census/02_genotyping/concat_vcfs/autosomes_output_filt_mindepth7_minqual30_AN22_chrrename.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out 11gorillas_unfilt_out

touch /global/scratch/users/rdekayne/gorilla_census/03_PCA/11_pca.done
```
Now run - `sbatch 03.1_PCA_unfilt.sh`
We get: 3823659 variants and 14 people pass filters and QC.

And then download the data from the cluster `scp rdekayne@hpc.brc.berkeley.edu:/global/scratch/users/rdekayne/gorilla_census/03_PCA/*.eigen* .`

Run the R script: `11_gorillas_pca.R`
```
#11 gorillas PCAs - De-Kayne 2024

#load the necessary packages
library(ggplot2)
library("ggrepel")
library(tidyverse)

setwd("/Users/rdk_mac_work/Dropbox/RishiMac3Work/gorilla_census/pca/")
background <- read.csv("../background/background_gorillas.csv", header = T, sep = ",")
str(background)
background <- background[,1:17]

eigenvec <- as.character("11gorillas_unfilt_out.eigenvec")
eigenval <- as.character("11gorillas_unfilt_out.eigenval")

pca <- read_table(eigenvec, col_names = FALSE)
pca <- pca[,-1]

#and load in the eigenval values
eigenval <- scan(eigenval)

names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

colour_plot <- rep(NA, length(pca$ind))
shape_plot <- rep(NA, length(pca$ind))

for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$namename) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$colour)
  shape_plot[i] <- as.character(col_name$shape)
}

pca <- as.data.frame(data.frame(pca, colour_plot, shape_plot))

pve <- data.frame(PC = 1:11, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

#cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("./11_gorilla_PC1PC2.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     cex = 1.75,
     main = "PC1 vs. PC2")
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)

pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 

dev.off()

tiff("./11_gorilla_PC1PC3.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), 
     cex = 1.75,
     main = "PC1 vs. PC3")
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)

pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 

dev.off()

tiff("./11_gorilla_PC2PC3.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), 
     cex = 1.75,
     main = "PC2 vs. PC3")
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 
dev.off()

###and again but this time just for the eigen files from the 2000 focal snps

eigenvec <- as.character("11gorillas_1000snp_panel_unfilt_out.eigenvec")
eigenval <- as.character("11gorillas_1000snp_panel_unfilt_out.eigenval")

pca <- read_table(eigenvec, col_names = FALSE)
pca <- pca[,-1]

#and load in the eigenval values
eigenval <- scan(eigenval)

names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

colour_plot <- rep(NA, length(pca$ind))
shape_plot <- rep(NA, length(pca$ind))

for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$namename) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$colour)
  shape_plot[i] <- as.character(col_name$shape)
}

pca <- as.data.frame(data.frame(pca, colour_plot, shape_plot))

pve <- data.frame(PC = 1:11, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

#cumulative sum of the percentage variance explained
cumsum(pve$pve)
tiff("./11_gorilla_2000_PC1PC2.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     cex = 1.75,
     main = "PC1 vs. PC2")
text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)

pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 
dev.off()

tiff("./11_gorilla_2000_PC1PC3.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), 
     cex = 1.75,
     main = "PC1 vs. PC3")
text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 

dev.off()


tiff("./11_gorilla_2000_PC2PC3.tiff", height=8, width=8, units="in", res=300)
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch = as.integer(pca$shape_plot),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), 
     cex = 1.75,
     main = "PC2 vs. PC3")
text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)

pop_list <- c("Mountain", "Eastern Lowland", "DRC", "Rwanda", "Uganda")
pop_list_shape <- c(15, 15, 16, 17, 18) 
pop_list_col <- c("darkslategray4", "goldenrod", "black", "black", "black")
legend('topright', 
       legend = pop_list, 
       col = pop_list_col,
       pch = as.integer(pop_list_shape),
       pt.cex = 1,
       cex = 1) 
dev.off()
```




