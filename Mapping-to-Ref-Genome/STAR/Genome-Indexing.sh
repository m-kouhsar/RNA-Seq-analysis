#!/bin/bash
#SBATCH -A  Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

######################################################################################################################

OutDir="/lustre/projects/Research_Project-T111004/Morteza/GenomeIndex/STAR/GRCh38.p14.PRI.GENECODE.R48"
Genome_fasta="/lustre/projects/Research_Project-T111004/Morteza/GenomeRef/GRCh38.p14.PRI.GENECODE.R48.fa"
Genome_annotation="/lustre/projects/Research_Project-T111004/Morteza/GenomeRef/GRCh38.p14.PRI.GENECODE.R48.gff3"

######################################################################################################################

module load STAR

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir $OutDir \
     --genomeFastaFiles $Genome_fasta \
     --sjdbGTFfile $Genome_annotation \
     --sjdbOverhang 150\
