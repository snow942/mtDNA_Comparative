#!/bin/bash -l

#SBATCH --job-name=SnowFreebayes

#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --exclusive=user

#SBATCH --output=./SlurmOutput/SnowFreebayes.%j.out  # STDOUT output file
#SBATCH --error=./SlurmOutput/SnowFreebayes.%j.err   # STDERR output file

#Load module(s)
enable_lmod

module load container_env
module load freebayes/1.3.1
module load vcftools/0.1.14
module load samtools/1.12
module load bcftools/1.10.2

#Documentation of version
freebayes --version
vcftools --version

#Downstream variant filter, variant must be present in 90% of community samples
cd /scratch/rsnow/mtDNA

BasCount=$(ls -d *Bas*.fq | wc -l)
ABasCount=$(ls -d *ABas*.fq | wc -l)
CBasCount=$(ls -d *CBas*.fq | wc -l)

BasFilter=$(echo $BasCount*.9 | bc)
ABasFilter=$(echo $ABasCount*.9 | bc)
CBasFilter=$(echo $CBasCount*.9 | bc)

#Running Freebayes for samples
#Samples manually moved to Sfa
cd /scratch/rsnow/Sfa

for bamfile in Sfa_*.bam
do
	echo "$bamfile"

	freebayes -f /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta -p 1 "$bamfile" > "$bamfile".vcf

#	if [[ "$bamfile" == *ABas* ]]; then
#		vcftools --vcf "$bamfile".vcf --minDP "$ABasFilter" --recode --minQ 20 --out "$bamfile"_filter
#
#	elif [[ "$bamfile" == *CBas* ]]; then
#		vcftools --vcf "$bamfile".vcf --minDP "$CBasFilter" --recode --minQ 20 --out "$bamfile"_filter
#
#	else
#		vcftools --vcf Sfa_ABas.merged.bam.vcf --minDP "$BasFilter" --recode --minQ 20 --out "$bamfile"_filter
#	fi
done

#Filter variants for a quality score greater than 20 and a minimum read depth of 90%
vcftools --vcf Sfa_ABas.merged.bam.vcf --minDP "$ABasFilter" --recode --minQ 20 --out Sfa_ABas_filter

vcftools --vcf Sfa_CBas.merged.bam.vcf --minDP "$CBasFilter" --recode --minQ 20 --out Sfa_CBas_filter

vcftools --vcf Sfa_Bas.merged.bam.vcf --minDP "$BasFilter" --recode --minQ 20 --out Sfa_Bas_filter


exit
