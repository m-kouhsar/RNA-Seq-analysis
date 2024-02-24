#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

######################################################################################################################

OutDir=./genome_index
Genome_fasta=./references/GRCh37.p13.genome.fa
Genome_annotation=./references/gencode.v19.annotation.gff3

######################################################################################################################

module load STAR

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir $OutDir \
     --genomeFastaFiles $Genome_fasta \
     --sjdbGTFfile $Genome_annotation \
     --sjdbOverhang 150\
