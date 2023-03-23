Mitogenome Comparative Analysis - ABas & CBas
========================================
Ryan Snow

Copmarative analysis of *Salarias fasciatas* mitochondrial assemblies from samples collected during the USS Albatross fieldwork (approximately 100 years ago) and current, contemporary assemblies. The intent of this study is to potentially identify any mitochondrial variants that may play a role in stress acclimation, associated with climate change responses. Initial thoughts behind this study were derived from several components; 
	1. mitochondrial plasticity and likelihood to adapt while mainting core functions 
	2. energy metabolism genes' role as an important factor in stress acclimation and adaptation responses across phyla 
	3. importance of oxygen availability in a warming ocean and the association with mitochondria. 

The data collected data was sequenced using a lcwgs technique. The low-coverage aspect of this study is regarded as a non-issue due to the abundance of mitochondria, and mitochondrial DNA, in cells, allowing for in-depth assemblies of the ~16,500 bp (mito)genome.

### ABas - USS Albatross Samples
### CBas - Contemporary Samples
### Bas - Combination of ABas and CBas Samples

# Steps for mitogenmoic comparative analysis:

1. [MUMmerQuery.sh](MUMmerQuery)
		
1. [Nucmer_mtDNA.sh](Nucmer_mtDNA)
	* The *.m files are to be used
		
1. delta2vcf
	* To be included in Nucmer_mtDNA bashscript
		
1. [SNPs2Counts.R](SNPs2Counts.R)
		
1. [mtDNAFilter.R](mtDNAFilter.R)
		
1. (SnpEff & SnpSift)[]
		
1. [Plink.sh](Plink)
		
1. [DiversityStats.R](DiversityStats.R)


# Steps for phylogenetic verification

1. MEGA...MUSCLE...


# Versions

* Nucmer - 
* R Studio - 
