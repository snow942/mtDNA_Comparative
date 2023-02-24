#Ryan Snow
#mtDNA_Comparative_Analysis

#Multiple sequence alignment completed via Nucmer, using the MUMmer package/software

library (dplyr)
library (devtools)
library (magrittr)
library (GenomicRanges)
library (knitr)
library (ggplot2)
library (tidyr)
library (rmarkdown)
library (rutilstimflutre)
library (ggbio)

#Read file(s)
Sfa_ABas.snps <- read.table (file="/Users/snow4/Desktop/MUMmer/Sfa_ABas_Query.fasta.snps")

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


#Determines number of total samples to implement a filter based on 90% of sample population
SampleCounts.df <- data.frame (table (Sfa_ABas.snps$SampleName))
TotalSamples <- nrow (SampleCounts.df)

Filter <- TotalSamples*0.9

#Counts the number of times a position appears in the dataset
Counts.df <- data.frame (table (Sfa_ABas.snps$P1))

#Filters Counts.df to show only variant positions in 90% of sample population
Pos.filtered <- as_tibble (
  Counts.df$Var1[(which (
    Counts.df$Freq >= Filter))])

VarPos <- data.frame (Pos.filtered$value)

Test <- GRanges (seqnames = Rle (VarPos$Pos.filtered.value),
                 ranges = IRanges (VarPos$Pos.filtered.value))

                 
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
  circle (Test, geom = "rect",
         aes (color = "red", fill = "red")) +
  circle (GRangeRef, geom = "ideo",
          aes (color = "gray", fill = "gray")) +
  # circle (GRangeRef, geom = "scale",
  #        scale.type = "sci", size = 1.5) +
  circle(Test, geom = "text",
         aes(label = seqnames), size = 1) 

Test
CirclePlot
