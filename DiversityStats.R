#Ryan Snow
#Diversity statistics

setwd ("C:/Users/snow4/Desktop/Mummerv4")

library (tidyverse)
library (dplyr)
library (ggplot2)
library (gtable)

####
#PCA
#https://speciationgenomics.github.io/pca/

Bas.PCA <- read_table ("Bas.Anno.plink.eigenvec", col_names = FALSE)
Bas.eigenval <- scan ("Bas.Anno.plink.eigenval")

Bas.PCA <- Bas.PCA [, -1]

names (Bas.PCA)[1] <- "Ind"
names (Bas.PCA)[2:ncol (Bas.PCA)] <- paste0 ("PC", 1:(ncol (Bas.PCA)-1))
Samples <- rep (NA, length (Bas.PCA$Ind))
Samples[grep ("Sfa-ABas", Bas.PCA$Ind)] <- "SfaABas"
Samples[grep ("Sfa-CBas", Bas.PCA$Ind)] <- "SfaCBas"
Samples[grep ("NC", Bas.PCA$Ind)] <- "Ref"

Bas.PCA <- as_tibble (data.frame (Bas.PCA, Samples)) 

Bas.PVE <- data.frame (PC = 1:20, Bas.PVE = Bas.eigenval/sum(Bas.eigenval)*100)

Bas.Hist <- ggplot (Bas.PVE, aes (PC, Bas.PVE)) + geom_bar (stat = "identity") 
Bas.Hist + ylab ("Percentage variance explained") + theme_classic ()

cumsum (Bas.PVE$Bas.PVE)

Bas.PCA.plot <- ggplot (Bas.PCA, aes (PC1, PC2, col = Samples)) + geom_point (size = 3) +
  scale_color_manual (values = c ("black", "green", "pink")) + 
  coord_equal () + theme_classic() +
  xlab (paste0("PC1 (", signif(Bas.PVE$Bas.PVE[1], 3), "%)")) + ylab(paste0("PC2 (", signif(Bas.PVE$Bas.PVE[2], 3), "%)"))
Bas.PCA.plot

####

library (vcfR)

ABas.Anno.vcf <- read.vcfR ("ABas.Anno.vcf")
CBas.Anno.vcf <- read.vcfR ("CBas.Anno.vcf")
Bas.Anno.vcf <- read.vcfR ("Bas.Anno.vcf")

#https://popgen.nescent.org/2015-12-15-microsatellite-differentiation.html#introduction

library (adegenet)
library (poppr)
library (hierfstat)
library (pegas)
library (mmod)
library (reshape2)

Bas.genind <- vcfR2genind(x = Bas.Anno.vcf)
Bas.Pop <- read.csv ("Bas_names.csv")
strata (Bas.genind) <- Bas.Pop
setPop (Bas.genind) <- ~Pop

Bas.PW <- pairwise_genetic_diff(vcf = Bas.SNPs2Counts, pops = Bas.genind@pop, method = "nei")
Bas.GD <- genetic_diff(vcf = Bas.SNPs2Counts, pops = Bas.genind@pop, method = "nei")
Bas.Fst <- pairwise.neifst(Bas.genind)
Bas.Stats <- basic.stats (Bas.genind, diploid = FALSE)
Bas.Hs <- Hs (Bas.genind)
Bas.Ho <- Ho (Bas.genind)

Bas.diff <- diff_stats (Bas.genind)
per.locus <- melt (Bas.diff$per.locus, varnames = c ("Locus", "Statistic"))
stats <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
Bas.glob <- data.frame (Statistic = stats, value = Bas.diff$global)

#AMOVA

Bas.dist <- dist (Bas.genind)
Bas.stra <- strata (Bas.genind)
Bas.AMOVA <- pegas::amova (Bas.dist ~ Pop, data = Bas.stra)

Bas.AMOVA
