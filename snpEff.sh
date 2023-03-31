#!/bin/bash -l

#SBATCH --job-name=MUM

#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --exclusive=user

#SBATCH --output=/home/rsnow/BashScripts/SlurmOutput/MUM.%j.err  # STDOUT output file
#SBATCH --error=/home/rsnow/BashScripts/SlurmOutput/MUM.%j.err   # STDERR output file

#Load Modules
enable_lmod

module load anaconda

conda deactivate
conda deactivate

conda activate Mummerv4

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/Mummerv4/

snpEff eff -c /home/rsnow/snpEff/snpEff.config Sfa_mtDNA ABas.SNPs2Counts.vcf > ABas.Anno.vcf

snpEff eff -c /home/rsnow/snpEff/snpEff.config Sfa_mtDNA CBas.SNPs2Counts.vcf > CBas.Anno.vcf

snpEff eff -c /home/rsnow/snpEff/snpEff.config Sfa_mtDNA Bas.SNPs2Counts.vcf > Bas.Anno.vcf

exit
