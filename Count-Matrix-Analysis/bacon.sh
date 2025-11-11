#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
####################################################################################################################

sumstat_file="./Results/NIH.snoRNA.limma.SV1.AD.C.csv"
effectsizes_col="logFC"
standarderrors_col="SE"
OutPrefix="./Results/NIH.limma.SV1"

ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Count-Matrix-Analysis"
#####################################################################################################################

Rscript ${ScriptDir}/bacon.R $sumstat_file $effectsizes_col $standarderrors_col $OutPrefix 

echo "All done!"


