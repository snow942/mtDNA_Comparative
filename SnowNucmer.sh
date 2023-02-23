#!/bin/bash -l

#SBATCH --job-name=MUM

#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --exclusive=user

#SBATCH --output=./SlurmOutput/SnowMUM.%j.out  # STDOUT output file
#SBATCH --error=./SlurmOutput/SnowMUM.%j.err   # STDERR output file

#Load Modules
enable_lmod

module load samtools/1.12

#Documentation of version
/cm/shared/apps/mummer/3.23/nucmer --version

#Creation of a list of the sample fastq to be used for while-loop
cd /scratch/rsnow/MUMmer


for mfafile in Sfa*
do

	/cm/shared/apps/mummer/3.23/nucmer --prefix="$mfafile" /home/rsnow/mtDNARef/Sfa_mtDNA_Ref.fasta "$mfafile"

#	/cm/shared/apps/mummer/3.23/mummerplot --prefix="$mfafile" --coverage --SNP "$mfafile".delta
#mummerplot commented out due to error in running

	/cm/shared/apps/mummer/3.23/show-coords -c "$mfafile".delta > "$mfafile".stdout

	/cm/shared/apps/mummer/3.23/show-snps -HlrT "$mfafile".delta > "$mfafile".snps

#	/cm/shared/apps/mummer/3.23/delta-filter -q "$mfafile".delta > "$mfafile".q.filter
#	/cm/shared/apps/mummer/3.23/show-snps -lr "$mfafile".?.filter > "$mfafile".filter.snps
done


exit
