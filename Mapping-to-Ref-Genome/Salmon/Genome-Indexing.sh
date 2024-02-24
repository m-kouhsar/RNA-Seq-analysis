#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

Genome_fasta=./Ref/GRCh38.primary_assembly.genome.fa.gz
Transcript_fasta=./Ref/gencode.v38.transcripts.fa.gz
OutDir=./Salmon/Genome_Index
Decoy_file=./Salmon/decoy.txt

module load Salmon

if [ ! -f $Decoy_file ]
then
	grep "^>" <(gunzip -c $Genome_fasta) | cut -d " " -f 1 > $Decoy_file
	sed -i.bak -e 's/>//g' $Decoy_file
fi

if [ ! -f ${OutDir}/genome_transcripts.fa.gz ]
then
	cat $Transcript_fasta $Genome_fasta > ${OutDir}/genome_transcripts.fa.gz
fi

salmon index -t ${OutDir}/genome_transcripts.fa.gz  \
	-i $OutDir \
	--decoys $Decoy_file \
  	-p 16 \
  	--gencode

echo "Indexing process has been done!"
