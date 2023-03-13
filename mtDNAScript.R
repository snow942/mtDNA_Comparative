#Ryan Snow
#mtDNA_Comparative_Analysis

library (dplyr)
library (devtools)
library (magrittr)
library (GenomicRanges)
library (knitr)
library (ggplot2)
library (tidyr)
library (rmarkdown)
library (tibble)


#Multiple sequence alignment completed via Nucmer, using the MUMmer package/software

#Read file(s)
Sfa_ABas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/SNPS/Sfa_ABas_Query.fasta.snps")
Sfa_CBas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/SNPS/Sfa_CBas_Query.fasta.snps")

Sfa_Bas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/SNPS/Sfa_Bas_Query.fasta.snps")

#Change column names to better represent data/incorporate metadata
#https://mummer.sourceforge.net/manual/#snps review output format for more info
colnames (Sfa_ABas.snps) <- c ("P1",
                               "SUB_Ref",
                               "SUB_Query",
                               "P2",
                               "BUFF",
                               "DIST",
                               "R",
                               "Q",
                               "LEN_R",
                               "LEN_Q",
                               "FRM1",
                               "FRM2",
                               "mtDNA_Rep",
                               "SampleName")

colnames (Sfa_CBas.snps) <- c ("P1",
                               "SUB_Ref",
                               "SUB_Query",
                               "P2",
                               "BUFF",
                               "DIST",
                               "R",
                               "Q",
                               "LEN_R",
                               "LEN_Q",
                               "FRM1",
                               "FRM2",
                               "mtDNA_Rep",
                               "SampleName")

colnames (Sfa_Bas.snps) <- c ("P1",
                              "SUB_Ref",
                              "SUB_Query",
                              "P2",
                              "BUFF",
                              "DIST",
                              "R",
                              "Q",
                              "LEN_R",
                              "LEN_Q",
                              "FRM1",
                              "FRM2",
                              "mtDNA_Rep",
                              "SampleName")

HighConfidenceVariants <- function (SnpsInput) {
#Used to call high confidence (90%) variants in .snps files
  ###SnpsInput <- Sfa_ABas.snps
  #Determines number of total samples to implement a filter based on 90% of sample population
  SampleCounts.df <- data.frame (table (SnpsInput$SampleName))
  TotalSamples <- nrow (SampleCounts.df)
  
  Filter <- TotalSamples*0.90
  
  #Counts the number of times a position appears in the dataset
  Positions.df <- data.frame (Var1=SnpsInput$P1, Ref=SnpsInput$SUB_Ref, Query=SnpsInput$SUB_Query)
  Counts.df <- data.frame (table (SnpsInput$P1))
  Positions.df$Var1 <- as.factor (Positions.df$Var1)
  
  To_Filter <- left_join (Counts.df, Positions.df, by=c("Var1"), multiple="first")
  colnames(To_Filter)[colnames(To_Filter) == "Var1"] ="POS"
  #Filters To_Filter to show only variant positions in 90% of sample population
  To_Filter.df <- data.frame (To_Filter)
  Output <- To_Filter.df %>% filter (Freq >= Filter)
    
  return (Output)
  
  }

#Calling high confidence variants in Sfa_ABas and Sfa_CBas files
Sfa_ABas.Filtered <- HighConfidenceVariants (Sfa_ABas.snps)
Sfa_CBas.Filtered <- HighConfidenceVariants (Sfa_CBas.snps)


print (Sfa_ABas.Filtered)
print (Sfa_CBas.Filtered)

ABas <- data.frame (Sfa_ABas.Filtered)
CBas <- data.frame (Sfa_CBas.Filtered)

#Compares ABas and CBas df to identify SNP positions unique to the population
Unique.ABas <- as_tibble (setdiff (ABas$Var1, CBas$Var1))
#Unique.CBas <- data.frame (setdiff (CBas$value, ABas$value))

            
####


#FakeDelta <- function (Filtered, Snps) {

#  Snps$P1 <- as.factor (Snps$P1)
#  Fake <- left_join (Filtered, Snps, by=c("Var1" = "P1"), multiple="first")

#  return (Fake)
  
#}

#ABas.Filtered.delta <- FakeDelta (ABas, Sfa_ABas.snps)
#CBas.Filtered.delta <- FakeDelta (CBas, Sfa_CBas.snps)

####

library ("vcfR")

ABas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_ABas_Query.fasta.vcf")
CBas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_CBas_Query.fasta.vcf")

ReadDepth <- function (vcf, Filtered.df) {
#Incorporate freq as QUAL scores for filtering of vcf
  
  vcf.df <- vcfR2tidy (vcf)
  Filtered.df$POS <- as.integer( as.character(Filtered.df$POS))
  #Changes class of df values to integer, where available
  
  merge.vcf <- full_join (vcf.df$fix, Filtered.df, by= "POS")
  #Merges the data from the vcf file and dataframe 
    #Uses "POS" column as like values
  merge.vcf <- merge.vcf %>% relocate (Freq, .after = "QUAL")
  #Moves Freq column position to after QUAL
  merge.vcf <- (merge.vcf[,!names(merge.vcf) %in%
                            c ("QUAL", "Ref", "Query", "ChromKey")])
  #Removes column with no signif value to exp
  colnames(merge.vcf)[colnames(merge.vcf) == "Freq"] ="QUAL"
  #Rename Freq column to QUAL
  merge.vcf$QUAL[is.na(merge.vcf$QUAL)] <- 0
  #All NAs in QUAL are changes to 0s
  INFO = c (NA)
  merge.vcf <- cbind(merge.vcf, INFO)
  #Appends INFO column that was lost due to unknown reasons
  vcf@fix <- as.matrix(merge.vcf)
  #Replaces old QUAL scores with new ones by replacing "fix" matrix in vcf
  
  return (vcf)
  }

ABas.DP.vcf <- ReadDepth (ABas.vcf, ABas)
CBas.DP.vcf <- ReadDepth (CBas.vcf, CBas)


write.vcf (ABas.DP.vcf, file = "ABas.Filter.vcf.gz", mask = FALSE, APPEND = FALSE)
write.vcf (CBas.DP.vcf, file = "CBas.Filter.vcf.gz", mask = FALSE, APPEND = FALSE)
#Files need to be written as zipped files