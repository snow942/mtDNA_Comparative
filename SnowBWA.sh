#!/bin/bash -l

#SBATCH --job-name=SnowBWA

#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --exclusive=user

#SBATCH --output=./SlurmOutput/SnowBWA.%j.out  # STDOUT output file
#SBATCH --error=./SlurmOutput/SnowBWA.%j.err   # STDERR output file

#Load Modules
enable_lmod

module load samtools/1.12
module load bwa/0.7.17

#Documentation of version
samtools --version
bwa mem --version

#Index reference for BWA
cd /scratch/rsnow/Sfa
bwa index /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta
samtools faidx /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta
touch Bas_FastqList.txt

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/mtDNA

for fastqname in Sfa-*.fq
do
	echo "$fastqname" >> Bas_FastqList.txt
done

#Running BWA on all files. Individual sample file vs Reference
cat Bas_FastqList.txt | while read name
do
	bwa mem /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta /scratch/rsnow/mtDNA/"$name" > "$name".sam

	samtools view -bt /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta "$name".sam > "$name".bam
	samtools fixmate -O bam "$name".bam "$name".fixmate.bam

	samtools sort -O bam -T "$name".sorted -o "$name".sorted.bam "$name".fixmate.bam

	samtools flagstat "$name".sorted.bam
done

#Merge files to respresent a community
samtools merge -f Sfa_Bas.merged.bam Sfa*.sorted.bam
samtools merge -f Sfa_CBas.merged.bam Sfa-CBas*.sorted.bam
samtools merge -f Sfa_ABas.merged.bam Sfa-ABas*.sorted.bam

samtools index -b Sfa_Bas.merged.bam



exit
