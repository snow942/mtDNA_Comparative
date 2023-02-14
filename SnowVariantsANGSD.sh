#!/bin/bash

#SBATCH --partition=main          # Partition (job queue)

#SBATCH --job-name=SnowVar       # Assign an short name to your job

#SBATCH --nodes=1                 # Number of nodes you require
#SBATCH --ntasks=1                # Total # of tasks across all nodes
#SBATCH --cpus-per-task=4         # Cores per task (>1 if multithread tasks)
#SBATCH --mem=200000                # Real memory (RAM) required (MB)
#SBATCH --time=48:00:00           # Total run time limit (HH:MM:SS)

#SBATCH --output=slurm.%N.%j.out  # STDOUT output file
#SBATCH --error=slurm.%N.%j.err   # STDERR output file (optional)

cd /scratch/rps109/BWA

module purge
module load intel/19.0.3 
module use /projects/community/modulefiles/
module load singularity/3.1.0
module load BWA/bwa-0.7.17-yc759.lua
module load samtools/1.8-gc563.lua


/home/rps109/angsd/angsd -bam bam.list -gl 1 -anc /home/rps109/mtDNA/Sfa_mtDNA_Ref.fasta -out ANGSD_Sfa -doSaf 1 -doMajorMinor 1 -doMaf 2 
#-SNP_pval 1e-6

###-gl 1 == SAMtools genotype likelihood model
###-anc == ancesteral state but can be supplied with ref
###-doSAF == Calculate the Site allele frequency assuming HWE
###-doMaf == Allele Frequency estimation



