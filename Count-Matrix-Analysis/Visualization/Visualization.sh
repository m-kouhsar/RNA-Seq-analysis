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
#########################################################################

results_file=
gene_col=args[2]
logFC_col=args[3]
pvalue_col=args[4]
pvalue_cut=args[5]
QQ_title=args[6]
Volcano_title=args[7]
OutPrefix=args[8]

ScriptDir=/lustre/projects/Research_Project-T111004/Morteza/github/RNA_Seq_analysis/Count-Matrix-Analysis/Visualization
