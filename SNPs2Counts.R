###Ryan Snow
###snps to vcf counts

library ("vcfR")

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
ABas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_ABas_Query.fasta.m.vcf")
ABas.names <- read.table ("ABas_names.txt")

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
CBas.vcf <- read.vcfR ("/Users/snow4/Desktop/Mummerv4/Sfa_CBas_Query.fasta.m.vcf")
CBas.names <- read.table ("CBas_names.txt")

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
Bas.names <- read.table ("Bas_names.txt")

#Merge snps to counts
SNPs2Counts <- function (SNPsFile, NamesFile, VcfFile) {
  
  mutable.vcf <- data.frame (VcfFile@fix)
  #Test
  #mutable.vcf <- data.frame (ABas.vcf@fix)
  #NamesFile <- ABas.names
  #SNPsFile <- Sfa_ABas.snps
  
  
  
  for (Name in NamesFile$V1) {
    mutable.vcf[Name] <- 0
  }

  for (location in NamesFile$V1) {
    for (number in mutable.vcf$POS){
      if (nrow (mutable.vcf[number == SNPsFile$P1 & location == SNPsFile$SampleName, ])) {
        print (number)
      
        rownumber <- which(mutable.vcf == number, arr.ind = TRUE)
        rownumber <- data.frame (rownumber)
      
        colnumber <- which(mutable.vcf == location, arr.ind = TRUE)
        colnumber <- data.frame (colnumber)
      
        locationnumber <- which(NamesFile == location, arr.ind = TRUE)
        locationnumber <- data.frame (locationnumber)
      
        mutable.vcf[rownumber$row, locationnumber$row+8] = (mutable.vcf[rownumber$row, locationnumber$row+8]) + 1
      }     
    }
  warnings ()
  }
  return (mutable.vcf)
}

ABas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_ABas.snps, NamesFile = ABas.names, VcfFile = ABas.vcf)
CBas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_CBas.snps, NamesFile = CBas.names, VcfFile = CBas.vcf)
Bas.SNPs2Counts <- SNPs2Counts(SNPsFile = Sfa_Bas.snps, NamesFile = Bas.names, VcfFile = Bas.vcf)

#Incorporating the new data into a vcf


#write.vcf (ABas.SNPs2Counts, file = "ABas.SNPs2Counts.vcf.gz", mask = FALSE, APPEND = FALSE)
#write.vcf (CBas.SNPs2Counts, file = "CBas.SNPs2Counts.vcf.gz", mask = FALSE, APPEND = FALSE)