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

conda activate Assemblytics

#
cd /scratch/rsnow/Mummerv4/

Assemblytics Sfa_ABas_Query.fasta.m.delta ABas.Assemblytics 1000 2 1000
Assemblytics Sfa_CBas_Query.fasta.m.delta CBas.Assemblytics 1000 2 1000
Assemblytics Sfa_Bas_Query.fasta.m.delta Bas.Assemblytics 1000 2 1000

exit
