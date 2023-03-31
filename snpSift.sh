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

SnpSift hwe Bas.Anno.vcf > Bas.Anno.hwe.vcf

SnpSift TsTv Bas.Anno.vcf > Bas.TsTv.txt

exit
