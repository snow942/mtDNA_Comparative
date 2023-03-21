#Ryan Snow
#mtDNA_Comparative_Analysis
####

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
Sfa_ABas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_ABas_Query.fasta.m.snps")
Sfa_CBas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_CBas_Query.fasta.m.snps")

Sfa_Bas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_Bas_Query.fasta.m.snps")

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

FilterVariants <- function (SnpsInput, FilterPercent) {
#Used to filter variants in .snps files based on read depth
  SnpsInput <- Sfa_ABas.snps
  #Determines number of total samples to implement a filter based on 90% of sample population
  SampleCounts.df <- data.frame (table (SnpsInput$SampleName))
  TotalSamples <- nrow (SampleCounts.df)
  
  #Filter
  Filter <- TotalSamples*(FilterPercent/100)
  
  #Counts the number of times a position appears in the dataset
  Positions.df <- data.frame (Var1=SnpsInput$P1, Ref=SnpsInput$SUB_Ref, Query=SnpsInput$SUB_Query)
  Counts.df <- data.frame (table (SnpsInput$P1))
  Positions.df$Var1 <- as.factor (Positions.df$Var1)
  
  To_Filter <- left_join (Counts.df, Positions.df, by=c("Var1"), multiple="first")
  colnames(To_Filter)[colnames(To_Filter) == "Var1"] ="POS"
  #Filters To_Filter to show only variant positions in 90% of sample population
  To_Filter.df <- data.frame (To_Filter)
  Output <- To_Filter.df %>% filter (Freq >= Filter)
  #Incorporates percentage of population with SNP per posiiton
  Output <- data.frame (Output, Percent = (as.integer(Output$Freq)/TotalSamples)*100)
    
  return (Output)
  
  }

#Calling high confidence variants in Sfa_ABas and Sfa_CBas files
ABas <- FilterVariants (Sfa_ABas.snps, 0)
CBas <- FilterVariants (Sfa_CBas.snps, 0)

#Compares filtered results of ABas and CBas populations
AC.compare <- full_join(ABas, CBas, by = "POS")

colnames (AC.compare) <- c ("Position", 
                            "ABas.Presence", "ABas.Ref", "ABas.Qry", "ABas.Percent",
                            "CBas.Presence", "CBas.Ref", "CBas.Qry", "CBas.Percent")
write.csv (AC.compare, file = "/Users/snow4/Desktop/Mummerv4/AC.compare.csv")

#Compares ABas and CBas df to identify SNP positions unique to the population
Unique.ABas <- as_tibble (setdiff (ABas$POS, CBas$POS))
Unique.CBas <- as_tibble (setdiff (CBas$POS, ABas$POS))

####

library ("vcfR")

ABas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_ABas_Query.fasta.1.vcf")
CBas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_CBas_Query.fasta.m.vcf")

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
  merge.vcf$QUAL[is.na(merge.vcf$QUAL)] <- 1
  #All NAs in QUAL are changes to 0s
  INFO = c (NA)
  #Appends INFO column that was lost due to unknown reasons
  merge.vcf <- cbind(merge.vcf, INFO)
  #Replaces old QUAL scores with new ones by replacing "fix" matrix in vcf
  vcf@fix <- as.matrix(merge.vcf)

  return (vcf)
  }

ABas.DP.vcf <- ReadDepth (ABas.vcf, ABas)
CBas.DP.vcf <- ReadDepth (CBas.vcf, CBas)

#Files need to be written as zipped files
write.vcf (ABas.DP.vcf, file = "/Users/snow4/Desktop/Mummerv4/ABas.m.Filter.vcf.gz", mask = FALSE, APPEND = FALSE)
write.vcf (CBas.DP.vcf, file = "/Users/snow4/Desktop/Mummerv4/CBas.m.Filter.vcf.gz", mask = FALSE, APPEND = FALSE)
