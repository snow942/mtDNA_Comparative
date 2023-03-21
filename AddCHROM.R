#Ryan Snow
#Add Chrom to mtDNA

library (vcfR)
library (dplyr)

ABas.noChrom <- read.vcfR ("/Users/snow4/Desktop/ClustalOmega/Aln.Sfa_ABas_Query.fasta.vcf")
CBas.noChrom <- read.vcfR ("/Users/snow4/Desktop/ClustalOmega/Aln.Sfa_CBas_Query.fasta.vcf")
Bas.noChrom <- read.vcfR ("/Users/snow4/Desktop/ClustalOmega/Aln.Sfa_Bas_Query.fasta.vcf")

ChromAdd <- function (vcf, Chrom = "") {
  #Incorporate Chrom value into vcf

  vcf <- ABas.noChrom
  Chrom = "NC_004412.1"
  vcf.df <- vcfR2tidy (vcf)
  SampleNum <- as.integer(dim(ABas.noChrom@gt)[1])
  Joining = 1
  Chrom_Num <- data.frame (REPLACEMENT = rep(Chrom, each = Joining), 
                           ChromKey = rep (Joining, each = Joining))
  
  merge.vcf <- left_join (x = vcf.df$fix, y = Chrom_Num, by = "ChromKey", multiple = "all")
  #Merges the data from the vcf file and Chrom
  
  merge.vcf <- merge.vcf %>% relocate (REPLACEMENT, .after = "CHROM")
  #Moves Freq column position to after QUAL
  merge.vcf <- (merge.vcf[,!names(merge.vcf) %in%
                            c ("CHROM", "ChromKey")])
  #Removes column with no signif value to exp
  colnames(merge.vcf)[colnames(merge.vcf) == "REPLACEMENT"] ="CHROM"
  #Rename REPLACEMENT column to CHROM
  INFO = c (NA)
  merge.vcf <- cbind(merge.vcf, INFO)
  #Appends INFO column that was lost due to unknown reasons
  vcf@fix <- as.matrix(merge.vcf)
  #Replaces old QUAL scores with new ones by replacing "fix" matrix in vcf
  
  return (vcf)
}

ABas.CHROM.vcf <- ChromAdd (ABas.noChrom, "NC_004412.1")
CBas.CHROM.vcf <- ChromAdd (CBas.noChrom, "NC_004412.1")
Bas.CHROM.vcf <- ChromAdd (Bas.noChrom, "NC_004412.1")

write.vcf (ABas.CHROM.vcf, file = "/Users/snow4/Desktop/ClustalOmega/ABas.CHROM.vcf.gz")
write.vcf (CBas.CHROM.vcf, file = "/Users/snow4/Desktop/ClustalOmega/CBas.CHROM.vcf.gz")
write.vcf (Bas.CHROM.vcf, file = "/Users/snow4/Desktop/ClustalOmega/Bas.CHROM.vcf.gz")
#Files need to be written as zipped files