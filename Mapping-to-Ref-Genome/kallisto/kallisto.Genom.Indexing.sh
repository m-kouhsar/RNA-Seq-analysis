#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=kallisto.GenomeIndex.%j.out

######################################################################################################################

OutPrefix=/lustre/projects/Research_Project-191391/Morteza/kallisto/Ref/GRCh38.p14
Genome_fasta="/lustre/projects/Research_Project-191391/Morteza/kallisto/Ref/gencode.v48.transcripts.fa"
Threads=15

######################################################################################################################

kallisto index -i "${OutPrefix}_kallistoIndex.idx" -t $Threads "$Genome_fasta"

echo "All done!"
