#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=15:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

######################################################################################################################

OutPrefix=/lustre/projects/Research_Project-191391/Morteza/kallisto/GRCh38.83.cdna
Genome_fasta=/lustre/projects/Research_Project-191391/Morteza/kallisto/Homo_sapiens.GRCh38.cdna.all.fa
Threads=15

######################################################################################################################

kallisto index -i "${OutPrefix}_kallistoIndex.idx" -t $Threads "$Genome_fasta"

