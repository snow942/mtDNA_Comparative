#!/bin/bash -l

#SBATCH --job-name=MUM

#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --exclusive=user

#SBATCH --output=/home/rsnow/BashScripts/SlurmOutput/MummerQuery.%j.out  # STDOUT output file
#SBATCH --error=/home/rsnow/BashScripts/SlurmOutput/MummerQuery.%j.err   # STDERR output file

#Load Modules

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/mtDNA

touch Sfa_ABas_Query.fasta
touch Sfa_CBas_Query.fasta
touch Sfa_Bas_Query.fasta

for fastaname1 in Sfa-ABas*.fa
do
        cat "$fastaname1" >> Sfa_ABas_Query.fasta
        sed -i "s/NC_004412.1/$fastaname1/g" Sfa_ABas_Query.fasta
done


for fastaname2 in Sfa-CBas*.fa
do
        cat "$fastaname2" >> Sfa_CBas_Query.fasta
        sed -i "s/NC_004412.1/$fastaname2/g" Sfa_CBas_Query.fasta
done


for fastaname3 in Sfa-*Bas*.fa
do
        cat "$fastaname3" >> Sfa_Bas_Query.fasta
        sed -i "s/NC_004412.1/$fastaname3/g" Sfa_Bas_Query.fasta
done


for seqname in *Query.fasta
do
        sed -i "s/_Ex1_L3_clmp.fp2_repr_fltrd.bam_mtseq.fa/ /g" "$seqname"
done

mv *Query.fasta /scratch/rsnow/Mummerv4/

exit
