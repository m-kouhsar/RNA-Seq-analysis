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
#########################################################################

counts_file=
pheno_file=
gFilter_min_count=5
gFilter_min_prop=0.75
model_main=~Phenotype+Age+Sex+RIN+PMI
model_null=~Age+Sex+RIN+PMI
var_interest=Phenotype
outliers=
OutPrefix=

ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Count-Matrix-Analysis/"
############################################################################

Rscripts "${ScriptDir}/SVA.R" "$counts_file" "$pheno_file" "$gFilter_min_count" "$gFilter_min_prop" "$model_main" "$model_null" "$var_interest" "$outliers" "$OutPrefix" 

echo "All done!"
