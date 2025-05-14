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
#SBATCH --output=Sleuth.%j.out

#################### Input Arguments ################################################################################

# SO_file: Sleuth Object file (created by 1.Sleuth.Read.R)
# lm_model: linear regression model to run the test
# factor_variables: Factor variable include condition variable in lm_model
# numeric_variables : Numerical varaibles in lm_model
# PCs: Number of principal components you want to add to the analysis as covariates
# OutPrefix: Results files/images prefix (can contains a directory)

#######################################################################################################################

SO_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/BDR.SO.rdat"
lm_model="~Phenotype+Age+Gender+RIN"
factor_variables="Phenotype,Gender"  
numeric_variables="Age,RIN"
PCs=0
OutPrefix="/lustre/projects/Research_Project-191391/Morteza/kallisto/BDR"
ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Mapping-to-Ref-Genome/kallisto/"

#######################################################################################################################

Rscript ${ScriptDir}/2.Sleuth.DEG.R $SO_file $lm_model $factor_variables $numeric_variables $PCs $OutPrefix 

