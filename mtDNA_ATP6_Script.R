#Ryan Snow
#mtDNA_Comparative_Analysis
#ATP6

library (dplyr)
library (devtools)
library (genbankr)
library (GenomicRanges)
#library (knitr)
#library (ggplot2)
library (tidyr)
library (rmarkdown)
#library (rutilstimflutre)
#library (ggbio)
#library (IRanges)

#Multiple sequence alignment completed via Nucmer, using the MUMmer package/software

#Read file(s)
Sfa_ABas.snps <- read.table (file="/Users/snow4/Desktop/ATP6/Sfa_ABas_Query.fasta.snps")
Sfa_CBas.snps <- read.table (file="/Users/snow4/Desktop/ATP6/Sfa_CBas_Query.fasta.snps")
Sfa_Bas.snps <- read.table (file="/Users/snow4/Desktop/ATP6/Sfa_Bas_Query.fasta.snps")

Ref_Anno <- readGenBank (file="/Users/snow4/Desktop/ATP6/ATP6.gb")

RefDb <- makeTxDbFromGenBank (Ref_Anno)

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
                               "FRM_1",
                               "FRM_2",
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
                               "FRM_1",
                               "FRM_2",
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
                               "FRM_1",
                               "FRM_2",
                               "mtDNA_Rep",
                               "SampleName")


HighConfidenceVariants <- function (SnpsInput) {
#Used to call high confidence (90%) variants in .snps files
  
  #Determines number of total samples to implement a filter based on 90% of sample population
  SampleCounts.df <- data.frame (table (SnpsInput$SampleName))
  TotalSamples <- nrow (SampleCounts.df)
  
  Filter <- TotalSamples*0.90
  
  #Counts the number of times a position appears in the dataset
  Counts.df <- data.frame (table (SnpsInput$P1))
  
  #Filters Counts.df to show only variant positions in 90% of sample population
  Pos.filtered <- as_tibble (
    Counts.df$Var1[(which (
      Counts.df$Freq >= Filter))])
  
  return (Pos.filtered)
  
  }

#Calling high confidence variants in Sfa_ABas and Sfa_CBas files
Sfa_ABas.Filtered <- HighConfidenceVariants (Sfa_ABas.snps)
Sfa_CBas.Filtered <- HighConfidenceVariants (Sfa_CBas.snps)


print (Sfa_ABas.Filtered)
print (Sfa_CBas.Filtered)

ABas <- data.frame (Sfa_ABas.Filtered)
CBas <- data.frame (Sfa_CBas.Filtered)

Unique.ABas <- data.frame (setdiff (ABas$value, CBas$value))
Unique.CBas <- data.frame (setdiff (CBas$value, ABas$value))


tx <- transcripts(RefDb, columns = "gene_id")



            
####

#Trying to create a circle plot that shows the location of each called var in the mitogenome

#https://taylorreiter.github.io/2019-05-11-Visualizing-NUCmer-Output/

FastaIndex <- read.table ("/Users/snow4/Desktop/mtDNARef/Sfa_mtDNA_Ref.fasta.fai", 
                          header = FALSE, stringsAsFactors = FALSE,
                          col.names = c ("Name", "Contig.Len", "Offset",
                                         "LineBases", "LineWidth"))

GRangeRef <- GRanges (seqnames = Rle (FastaIndex$Name),
                      ranges = IRanges (start = rep (1, nrow (FastaIndex)),
                                        end = FastaIndex$Contig.Len))

seqlengths (GRangeRef) <- FastaIndex$Contig.Len
genome (GRangeRef) <- "Reference"
GRangeRef <- sortSeqlevels (GRangeRef)
GRangeRef <- sort (GRangeRef)



CirclePlot <- ggbio () +
  circle (ABas, geom = "rect",
         aes (color = "red", fill = "red")) +
  circle (GRangeRef, geom = "ideo",
          aes (color = "gray", fill = "gray")) +
  # circle (GRangeRef, geom = "scale",
  #        scale.type = "sci", size = 1.5) +
  circle(ABas$value, geom = "text",
         aes(label = seqnames), size = 1) 

CirclePlot
