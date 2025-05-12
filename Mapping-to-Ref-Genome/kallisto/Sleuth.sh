#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out

#######################################################################################################################
kallisto_res_dir="/lustre/projects/Research_Project-191391/Morteza/kallisto/UKBBN/"
# pheno_file is a csv file which must contains the following columns:
#        sample: Samples ID
#        path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
#        All the variables that are included in lm_model must be represented by a column with the same name
pheno_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/UKBBN.Phenotype.Labeled.csv"
# Target mapping file contains information (eg. gene ID and type) about the transcripts
target_map_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/Human.Transcript.Ensemble.GRCh38.txt"
lm_model="~Phenotype+Age+Gender+RIN"
# Factor variable include condition variable in lm_model
factor_variables="Phenotype,Gender"  
# Numerical varaibles in lm_model
numeric_variables="Age,RIN"
# Number of Principal Components to add to the lm model (can be 0)
PCs=10
OutPrefix="/lustre/projects/Research_Project-191391/Morteza/kallisto/UKBBN"

ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Mapping-to-Ref-Genome/kallisto/"

#######################################################################################################################

Rscript ${ScriptDir}/Sleuth.R $kallisto_res_dir $pheno_file $target_map_file $lm_model $factor_variables $numeric_variables $PCs $OutPrefix $ScriptDir 

