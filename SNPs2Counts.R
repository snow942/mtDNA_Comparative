###Ryan Snow
###snps to vcf counts

library ("vcfR")
library ("dplyr")

#File import
setwd("C:/Users/snow4/Desktop/Mummerv4")

Sfa_ABas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_ABas_Query.fasta.m.snps")
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
ABas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/ABas.m.Filter.vcf")
ABas.names <- read.table ("Sfa_ABas_Query.fasta.txt")

Sfa_CBas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_CBas_Query.fasta.m.snps")
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
CBas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/CBas.m.Filter.vcf")
CBas.names <- read.table ("Sfa_CBas_Query.fasta.txt")

Sfa_Bas.snps <- read.table (file="/Users/snow4/Desktop/Mummerv4/Sfa_Bas_Query.fasta.m.snps")
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
Bas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_Bas_Query.fasta.m.vcf")
Bas.names <- read.table ("Sfa_Bas_Query.fasta.txt")

#Merge snps to counts
SNPs2Counts <- function (SNPsFile, NamesFile, VcfFile) {
  
  OutputVCF <- VcfFile
  mutable.vcf <- data.frame (VcfFile@fix)
  
  #Test
  #mutable.vcf <- data.frame (ABas.vcf@fix)
  #NamesFile <- ABas.names
  #SNPsFile <- Sfa_ABas.snps
  
  for (Name in NamesFile$V1) {
    mutable.vcf[Name] <- 0
    #Creates columns with sample name as header and 0 as data
  }
  
  for (location in NamesFile$V1) {
  #Cycles through each sample name
    
    for (number in mutable.vcf$POS){
    #Cycles through each SNP position called
      
      if (nrow (mutable.vcf[number == SNPsFile$P1 & location == SNPsFile$SampleName, ])) {
      #Verifies SNP position for each sample by comparing to SNPs file from Mummer
        #print (number)
      
        #Which row is the called SNP in the working file
        rownumber <- which(mutable.vcf == number, arr.ind = TRUE)
        rownumber <- data.frame (rownumber)
      
        #Which column is the called sample in the working file
        colnumber <- which(mutable.vcf == location, arr.ind = TRUE)
        colnumber <- data.frame (colnumber)
      
        #Which column is the called sample in the working file
        locationnumber <- which(NamesFile == location, arr.ind = TRUE)
        locationnumber <- data.frame (locationnumber)
      
        mutable.vcf[rownumber$row, locationnumber$row+8] = (mutable.vcf[rownumber$row, locationnumber$row+8]) + 1
        #Replaces 0 with 1 in SNP and Sample cell of interest (per cycle)
      }     
    }
  warnings ()
  #Reports any warnings
  }
  
 
  
  FORMAT = c ("GT")
  #Create new, updated FORMAT column
  merge.vcf <- cbind(mutable.vcf, FORMAT)
  #Appends FORMAT column to working dataframe
  merge.vcf <- merge.vcf %>% relocate (FORMAT)
  #Moves new FORMAT column position to first column position
  merge.vcf <- (merge.vcf[,!names(merge.vcf) %in%
                            c ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")])
  #Removes unwanted info and columns ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  OutputVCF@gt <- as.matrix(merge.vcf)
  #Replaces old QUAL scores with new ones by replacing "fix" matrix in vcf
  
  return (OutputVCF)
}



ABas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_ABas.snps, NamesFile = ABas.names, VcfFile = ABas.vcf)
CBas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_CBas.snps, NamesFile = CBas.names, VcfFile = CBas.vcf)
Bas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_Bas.snps, NamesFile = Bas.names, VcfFile = Bas.vcf)

#Incorporating the new data into a vcf
write.vcf (ABas.SNPs2Counts, file = "ABas.SNPs2Counts.vcf.gz", mask = FALSE, APPEND = FALSE)
write.vcf (CBas.SNPs2Counts, file = "CBas.SNPs2Counts.vcf.gz", mask = FALSE, APPEND = FALSE)
write.vcf (Bas.SNPs2Counts, file = "Bas.SNPs2Counts.vcf.gz", mask = FALSE, APPEND = FALSE)

