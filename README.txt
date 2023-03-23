mtDNA Comparative Analysis - ABas & CBas
========================================
Ryan Snow

Copmarative analysis of *Salarias fasciatas* mitochondrial assemblies from samples collected during the USS Albatross fieldwork (approximately 100 years ago) and current, contemporary assemblies. The intent of this study is to potentially identify any mitochondrial variants that may play a role in stress acclimation, associated with climate change responses. Initial thoughts behind this study were derived from several components; 1. mitochondrial plasticity and likelihood to adapt while mainting core functions 2. energy metabolism genes' role as an important factor in stress acclimation and adaptation responses across phyla 3. importance of oxygen availability in a warming ocean and the association with mitochondria. 

The data collected data was sequenced using a lcwgs technique. The low-coverage aspect of this study is regarded as a non-issue due to the abundance of mitochondria, and mitochondrial DNA, in cells, allowing for in-depth assemblies of the ~16,500 bp (mito)genome.

## ABas - USS Albatross Samples
## CBas - Contemporary Samples
## Bas - Combination of ABas and CBas Samples

# Steps for mitogenmoic comparative analysis:

1. (MUMmerQuery)[MUMmerQuery.sh]
		|
2a. (Nucmer_mtDNA)[Nucmer_mtDNA.sh]
	The *.m* files are to be used
		|
2b. delta2vcf
	To be included in Nucmer_mtDNA bashscript
		|
3. (SNPs2Counts.R)[SNPs2Counts.R]
		|
4. (mtDNAFilter.R)[mtDNAFilter.R]
		|
5. (SnpEff & SnpSift)[]
		|
6. (Plink)[Plink.sh]
		|
7. (DiversityStats.R)[DiversityStats.R]


# Steps for phylogenetic verification

1. MEGA...MUSCLE...


# Versions

* Nucmer - 
* R Studio - 
