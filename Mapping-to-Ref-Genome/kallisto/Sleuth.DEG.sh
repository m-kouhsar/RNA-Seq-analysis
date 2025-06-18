#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=2:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=Sleuth.DEG.%j.out

#################### Input Arguments ################################################################################

# SO_file: Sleuth Object file (created by 1.Sleuth.Read.R)
# lm_model: linear regression model to run the test
# PCs: Number of principal components you want to add to the analysis as covariates
# SVs: Number of surrogate variables you want to add to the analysis as covariates
# OutPrefix: Results files/images prefix (can contains a directory)

#######################################################################################################################

SO_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/Results/NIH.111.SO.rdat"
lm_model="~Phenotype+Age+Gender+RIN"
PCs=10
SVs=0
OutPrefix="/lustre/projects/Research_Project-191391/Morteza/kallisto/Results/NIH.PC10"


ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Mapping-to-Ref-Genome/kallisto/"

#######################################################################################################################

Rscript ${ScriptDir}/Sleuth.DEG.R $SO_file $lm_model $PCs $SVs $OutPrefix 

