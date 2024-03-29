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

#Documentation of version
plink --version

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/Mummerv4/


plink --vcf /scratch/rsnow/Mummerv4/Bas.Anno.hwe.vcf --make-bed -chr-set 1 no-xy --allow-extra-chr --out Bas_plink

plink --bfile Bas_plink --pca --chr-set 1 no-xy --allow-extra-chr --out Bas_plink

exit
