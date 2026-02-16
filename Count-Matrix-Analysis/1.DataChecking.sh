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

counts_file="./Raw/NIH/NIH_miRNA.mirdeep2.count.tsv"
pheno_file="./NIH/NIH.Phenotype.csv"
trait_var="Phenotype"
factor_vars="Sex,RINcat,PMICat,Plate"
numeric_vars="Age"
normalize_method="cpm"
lib_size_threshold="1e+6"
gFilter_min_count=5
gFilter_min_prop=0.75
OutPrefix="Results/NIH"

ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Count-Matrix-Analysis"
##########################################################################

Rscript $ScriptDir/DataChecking.R $ScriptDir $counts_file $pheno_file $trait_var $factor_vars $numeric_vars $normalize_method $lib_size_threshold $gFilter_min_count $gFilter_min_prop $OutPrefix

echo "All done!"