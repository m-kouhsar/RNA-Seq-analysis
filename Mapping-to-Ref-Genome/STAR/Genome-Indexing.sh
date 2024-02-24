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

OutDir=/lustre/projects/Research_Project-191391/STAR/references/genome_index37
Genome_fasta=/lustre/projects/Research_Project-191391/STAR/references/genome_fasta/GRCh37.p13.genome.fa
Genome_annotation=/lustre/projects/Research_Project-191391/STAR/references/genome_anno/gencode.v19.annotation.gff3

######################################################################################################################

module load STAR

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir $OutDir \
     --genomeFastaFiles $Genome_fasta \
     --sjdbGTFfile $Genome_annotation \
     --sjdbOverhang 150\
