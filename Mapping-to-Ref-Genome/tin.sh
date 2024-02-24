#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=250:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

############################################################
##    https://rseqc.sourceforge.net/#tin-py 
############################################################
BamFileDir=./BamFiles
RefGenomeBed=./Reference/gencodev38.bed

module load RSeQC

tin.py -i $BamFileDir -r $RefGenomeBed

