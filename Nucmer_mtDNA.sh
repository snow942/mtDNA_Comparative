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

conda activate Mummer4

#Documentation of version
nucmer --version

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/Mummerv4


for mfafile in Sfa*
do

	nucmer --prefix="$mfafile" /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta "$mfafile"

	delta-filter -1 "$mfafile".delta | delta2vcf > "$mfafile".1.vcf
	delta-filter -m "$mfafile".delta | delta2vcf > "$mfafile".m.vcf

	delta-filter -1 "$mfafile".delta > "$mfafile".1.delta
	delta-filter -m "$mfafile".delta > "$mfafile".m.delta

	show-snps -HlrT "$mfafile".1.delta > "$mfafile".1.snps
	show-snps -HlrT "$mfafile".m.delta > "$mfafile".m.snps

	show-diff -q "$mfafile".1.delta > "$mfafile".1.report
	show-diff -q "$mfafile".m.delta > "$mfafile".m.report

done

exit
