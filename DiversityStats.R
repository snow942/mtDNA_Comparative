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

Bas.PCA <- read_table ("Bas_plink.eigenvec", col_names = FALSE)
Bas.eigenval <- scan ("Bas_plink.eigenval")

Bas.PCA$X2 <- str_c (Bas.PCA$X1, Bas.PCA$X2)
Bas.PCA <- Bas.PCA [, -1]

names (Bas.PCA)[1] <- "Ind"
names (Bas.PCA)[2:ncol (Bas.PCA)] <- paste0 ("PC", 1:(ncol (Bas.PCA)-1))
Samples <- rep (NA, length (Bas.PCA$Ind))
Samples[grep ("Sfa-ABas", Bas.PCA$Ind)] <- "SfaABas"
Samples[grep ("Sfa-CBas", Bas.PCA$Ind)] <- "SfaCBas"
#Samples[grep ("NC", Bas.PCA$Ind)] <- "Ref"

Bas.PCA <- as_tibble (data.frame (Bas.PCA, Samples)) 

Bas.PVE <- data.frame (PC = 1:20, Bas.PVE = Bas.eigenval/sum(Bas.eigenval)*100)

Bas.Hist <- ggplot (Bas.PVE, aes (PC, Bas.PVE)) + geom_bar (stat = "identity") 
Bas.Hist + ylab ("Percentage variance explained") + theme_classic ()

cumsum (Bas.PVE$Bas.PVE)

Bas.PCA.plot <- ggplot (Bas.PCA, aes (PC1, PC2, col = Samples)) + 
  geom_point (aes (shape = Samples), size = 3) +
  #scale_color_manual (values = c ("green", "pink")) + 
  coord_equal () + theme_classic() +
  xlab (paste0("PC1 (", signif(Bas.PVE$Bas.PVE[1], 3), "%)")) + ylab(paste0("PC2 (", signif(Bas.PVE$Bas.PVE[2], 3), "%)"))

Bas.PCA.plot

####

library (vcfR)

#ABas.Anno.vcf <- read.vcfR ("ABas.Anno.vcf")
#CBas.Anno.vcf <- read.vcfR ("CBas.Anno.vcf")
Bas.Anno.vcf <- read.vcfR ("Bas.Anno.hwe.vcf")

#https://popgen.nescent.org/2015-12-15-microsatellite-differentiation.html#introduction

library (adegenet)
library (poppr)
library (hierfstat)
library (pegas)
library (mmod)
library (reshape2)

Bas.genind <- vcfR2genind (x = Bas.Anno.vcf)

Bas.Pop <- as_tibble (x=Bas.PCA$Ind)
Bas.Pop <- Bas.Pop %>% add_column (Bas.PCA$Samples)
colnames (Bas.Pop) <- c ("Indv", "Pop")

strata (Bas.genind) <- Bas.Pop
setPop (Bas.genind) <- ~Pop

#Bas.PW <- pairwise_genetic_diff(vcf = Bas.Anno.vcf, pops = Bas.genind@pop, method = "nei")
#Bas.PW
#Bas.GD <- genetic_diff(vcf = Bas.Anno.vcf, pops = Bas.genind@pop, method = "nei")
#Bas.GD

Bas.Fst <- pairwise.WCfst(Bas.genind)
Bas.Fst

Bas.Stats <- basic.stats (Bas.genind, diploid = FALSE)
Bas.Stats
#Bas.Hs <- Hs (Bas.genind)
#Bas.Hs

#Bas.Ho <- Ho (Bas.genind)
#Bas.Ho

#Bas.diff <- diff_stats (Bas.genind)
#Bas.diff
#per.locus <- melt (Bas.diff$per.locus, varnames = c ("Locus", "Statistic"))
#stats <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
#Bas.glob <- data.frame (Statistic = stats, value = Bas.diff$global)
#Bas.glob

#Principal component(s) matched with assoc. sample ID in a tibble
#Created for Fst
PCA_Col_Comb <- function (Sample_Col_Input, PCA_Col_input, header="PC", max_lim, min_lim) {
  
  Pop.PCA <- as_tibble (x=Sample_Col_Input)
  PCA_Col_Out <- as_tibble ("")
  #Creation of empty tibble
  
  for (x in PCA_Col_input) {
    if (x >= as.numeric (max_lim)) {
      PCA_Col_Out [nrow (PCA_Col_Out) + 1,] <- c ("Group1")
    }
    #If sample PC is greater than defined upper limit place into group1
    if (x <= as.numeric (min_lim)) {
      PCA_Col_Out [nrow (PCA_Col_Out) + 1,] <- c ("Group2")
    }
    #If sample PC is less than defined lower limit place into group2
    else if (x > as.numeric (min_lim) & x < as.numeric (max_lim)) {
      PCA_Col_Out [nrow (PCA_Col_Out) + 1,] <- c ("Group3")
    }
    #If sample PC is bewteen defined limits place into group3
  }
  print (min_lim)
  PCA_Col_Out = PCA_Col_Out[-1,]
  Pop.PCA <- Pop.PCA %>% add_column (PCA_Col_Out$value)
  colnames (Pop.PCA) <- c ("Indv", header)
  
  return (Pop.PCA)
}

#PC1
Bas.Pop.PCA1 <- PCA_Col_Comb (Sample_Col_Input = Bas.PCA$Ind, PCA_Col_input = Bas.PCA$PC1, 
                             header = "PC1", max_lim = 0.3, min_lim = 0.1)

strata (Bas.genind) <- Bas.Pop.PCA1
setPop (Bas.genind) <- ~PC1
Bas.PC1.Fst <- pairwise.WCfst(Bas.genind)
Bas.PC1.Fst

#PC2
Bas.Pop.PCA2 <- PCA_Col_Comb (Sample_Col_Input = Bas.PCA$Ind, PCA_Col_input = Bas.PCA$PC2, 
                             header = "PC2", max_lim = 0.2, min_lim = (-0.1))

strata (Bas.genind) <- Bas.Pop.PCA2
setPop (Bas.genind) <- ~PC2
Bas.PC2.Fst <- pairwise.WCfst(Bas.genind)
Bas.PC2.Fst

#AMOVA

Bas.dist <- dist (Bas.genind)
Bas.stra <- strata (Bas.genind)
Bas.AMOVA <- pegas::amova (Bas.dist ~ Pop, data = Bas.stra)

Bas.AMOVA

#Tajima

Bas.DNABin <- vcfR2DNAbin (x = Bas.Anno.vcf)
Bas.Taj <- tajima.test(x = Bas.DNABin)
Bas.Taj
