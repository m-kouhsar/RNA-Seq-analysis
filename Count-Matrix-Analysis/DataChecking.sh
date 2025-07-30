#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=01:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#########################################################################

counts_file="./Raw/Project_11556.STAR.count.reverse.tsv"
pheno_file="./Raw/Pheno_Ascribed.csv"
trait_vars="Trait"
factor_vars="Sex,Plate"
numeric_vars="Age,RIN"
OutPrefix="./Results/Ascribe_Study"

ScriptDir="./RNA_Seq_analysis/Count-Matrix-Analysis/"
##########################################################################

Rscript $ScriptDir/DataChecking.R $counts_file $pheno_file $trait_vars $factor_vars $numeric_vars $OutPrefix

echo "All done!"